/*
   Example compound:
     a-b=a-c-a-a=b-c 
     
   Reactions center around bonds; the
   reaction between two molecules happens
   at the bonds with the most similar
   chemical environments
   
   Decomposition is also possible, and tends to create a double bond
   
   There are two general motifs:
  
   -a=b- + -c-d-  -> -a-d- + -c-b-
   -a-b- <-> =a + b-
*/

var maxLength = 50;
var OFFDIAG = 0.1; // How much transformation on chemical identity is performed when you traverse a bond
var LOCALITY = 1.0; //0.5; // How quickly chemical environment falls off (looks like exp(-LOCALITY*x))
var TEMPERATURE = 0.1; // Reactions occur with rate ~ exp(-deltaE/TEMPERATURE)

var chemVecSize = 8; // This determines the richness of the 'chemical identity' vector. If this is less than the number of kinds of atoms, there is degeneracy/interchangeability

var compoundHash = []; /* This stores the energy and chemical environment for compounds we've previously encountered, to speed things up */

var chemBath = [];

var decompVector = new Vector(chemVecSize); // This particular vector must be within a certain distance of the local chemical identity for decomposition processes to occur at that bond

var bondEnergies = []; // Stores the energy of each type of single/double bond (A-A, B=C, etc)
var bondMatrix = [ new Matrix(chemVecSize), new Matrix(chemVecSize) ]; // This matrix multiplies the chemical identity of an atom type when you traverse a bond
var chemID = []; // This stores the chemical identity of each atom type

var rxnList = []; // This stores the most recent reactions for analysis

/* Initialize the bond matrices with random values */
function setupBondMatrix()
{
   for (var k=0;k<2;k++)
     for (var j=0;j<chemVecSize;j++)
       for (var i=0;i<chemVecSize;i++)
       {
          bondMatrix[k].x[i+j*chemVecSize] = 
	     LOCALITY*( (i==j) + OFFDIAG*(2*Math.random()-1));
       }
}

/* Initialize the chemical identities with random values */
function setupChemID()
{
   var charArray = ['A','B','C','D'];
   
   for (var j=0;j<4;j++)
   {
		chemID[charArray[j]] = new Vector(chemVecSize);

		for (var i=0;i<chemVecSize;i++)
		{
			chemID[charArray[j]].x[i]=2*Math.random()-1;
		}
	}
	
   for (var i=0;i<chemVecSize;i++)
      decompVector.x[i]=2*Math.random()-1;
      
	decompVector = normalize(decompVector);
}

/* Initialize the bond energies. Bond energies are symmetric in this 
 * model (but don't have to be in general), and double bonds are simply 
 * 'slightly deeper' than single bonds in a homogeneous way
 */
function setupBondEnergies()
{
   // A,B,C,D
   
   bondEnergies["A-A"]=-0.5;
   bondEnergies["A-B"]=-0.7;
   bondEnergies["A-C"]=-0.1;
   bondEnergies["A-D"]=-0.8;

   bondEnergies["B-B"]=-1.0;
   bondEnergies["B-C"]=-1.5;
   bondEnergies["B-D"]=-1.3;

   bondEnergies["C-C"]=-0.6;
   bondEnergies["C-D"]=-1.3;
   
   bondEnergies["D-D"]=1.0;

   var charArray = ['A','B','C','D'];
   
   for (var i=0;i<charArray.length;i++)
   {
      for (var j=i;j<charArray.length;j++)
      {
         var base=bondEnergies[charArray[i]+'-'+charArray[j]];
         bondEnergies[charArray[i]+'='+charArray[j]]=1.2*base;
         bondEnergies[charArray[j]+'='+charArray[i]]=1.2*base;
         bondEnergies[charArray[j]+'-'+charArray[i]]=base;	 
      }
   }
}

/* This defines the 'Chemical' object, which has a sequence, an energy,
 * and a chemical environment for each bond */
 
function Chemical(sequence)
{
	this.sequence = sequence; 

	// Check the hash to see if we have ever made this compound before
	var hashEntry = compoundHash[sequence];
		
	// If we have, use pre-computed values
	if (hashEntry) 
	{
		this.energy = hashEntry.energy;
		this.chemEnv = hashEntry.chemEnv;
	}
	else // Otherwise compute them and store them in the hash
	{
		this.energy = this.getEnergy();
   
		this.chemEnv = [];
	
		// Only have an environment around a bond
		for (var i=1;i<sequence.length-1;i+=2)
			this.chemEnv[i] = this.getChemEnv(i);
			
		compoundHash[sequence]={ energy: this.energy, chemEnv: this.chemEnv };
	}
}

/* This calculates the energy of a compound. In this model, it is purely
 * based on counting bonds. To generalize this, modify getEnergy
 */
Chemical.prototype.getEnergy = function()
{
	var energy=0;
   
	for (var i=1;i<this.sequence.length-1;i+=2)
	{
		energy+=bondEnergies[this.sequence.substr(i-1,3)];
	}
   
	return energy;
}

/* This calculates the chemical environment at a bond at position 'pos'
 * Note that 'pos' must be a bond location (e.g. an odd index) or this
 * will behave very poorly, as it assumes that every 2 letters away
 * from 'pos' will be a bond.
 */
Chemical.prototype.getChemEnv = function(pos)
{
   var chemEnv = new Vector(chemVecSize);
   var cMatrix = new Matrix(chemVecSize);
   
   for (var i=pos;i<this.sequence.length-1;i+=2)
   {
      var bond_id=0; if (this.sequence[i]=='=') bond_id=1;
      
      cMatrix=mulMatrix(bondMatrix[bond_id],cMatrix);
      chemEnv=vectAdd(chemEnv,matrixMulVector(cMatrix,chemID[this.sequence[i+1]]));
   }
   
   cMatrix.identity();
   
   for (var i=pos;i>=1;i-=2)
   {
      var idx=0;
      
      cMatrix=mulMatrix(bondMatrix[bond_id],cMatrix);
      chemEnv=vectAdd(chemEnv,matrixMulVector(cMatrix,chemID[this.sequence[i-1]]));
   }
   
   chemEnv = normalize(chemEnv);
   return chemEnv;
}

/* This function attempts to find a valid decomposition reaction and
 * evaluates whether it occurs via a Boltzmann exp(-deltaE/TEMPERATURE)
 * term. The 'best' site is always used (where 'best' is determined by
 * chemical environment of the bond, not the energy)
 *
 * The reaction is of the form: 
 *    -a-b- -> =a + b-
 * 
 * Note that this is direction-dependent (in this model). The following
 * reaction is not tested for:
 *    -a-b- -> -a + b=
 */

function tryDecompose(chem1)
{
   var outcome = {
      prod1: null,
      prod2: null
   };
   
   var bestdist=1e9;
   var site = null;

	// Cannot get lone atoms
   for (var i=3;i<chem1.sequence.length-3;i++)
   {
		if (chem1.sequence[i]=='-') // viable site
		{
			var dchem=vectAdd(chem1.chemEnv[i],scalMult(-1,decompVector));
			var dist = dot(dchem,dchem);
					
			// Only choose the one best reaction as far as chemical environment
			if ((dist<bestdist)&&(dist<2.0)) 
			{
				bestdist=dist;
				site = { si: i };
			}
		}
	}
	
	if (site)
	{
		// Example: a-b-c-d -> a=b + c-d ; site.si = 3
		outcome.prod1 = new Chemical(chem1.sequence.slice(0,site.si-2) +
									"=" + chem1.sequence[site.si-1]);
									
		outcome.prod2 = new Chemical(chem1.sequence.slice(site.si+1,chem1.sequence.length));

		var deltaE = (outcome.prod1.energy + outcome.prod2.energy) -
					 (chem1.energy);
	   
		if (Math.random()<Math.exp(-deltaE/TEMPERATURE))
			return outcome;
		else return null;
	}
	
	return null;
}

/* This is the reverse of the decomposition reaction
 * =a + b- -> -a-b-
 * 
 * It does however check for the opposite shape
 * -a + b= -> -a-b-
 */
function tryRecompose(chem1, chem2)
{
	var outcome = {
		prod1: null,
	};
	var site = null;
	var bestdist = 0
      
	if (chem1.sequence[chem1.sequence.length-2]=='=')
	{
		var dchem=vectAdd(chem1.chemEnv[chem1.sequence.length-2],
					scalMult(-1,chem2.chemEnv[1]));
		var dist = dot(dchem,dchem);
					
		// Only choose the one best reaction as far as chemical environment
		if (dist>bestdist) 
		{
			bestdist=dist;
			site = { si: chem1.sequence.length-2, sj: 1 };
		}
	}
	
	if (chem1.sequence[1]=='=')
	{
		var dchem=vectAdd(chem1.chemEnv[1],
					scalMult(-1,chem2.chemEnv[chem2.sequence.length-2]));
		var dist = dot(dchem,dchem);
					
		// Only choose the one best reaction as far as chemical environment
		if (dist>bestdist) 
		{
			bestdist=dist;
			site = { si: 1, sj: chem2.sequence.length-2 };
		}
	}
	
	if (site)
	{
		var newseq;
		
		if (site.si == 1)
			newseq = chem2.sequence + '-' + chem1.sequence[0] + 
					'-' + chem1.sequence.slice(2,chem1.sequence.length);
		else
			newseq = chem1.sequence.slice(0,chem1.sequence.length-2) + 
			'-' + chem1.sequence[chem1.sequence.length-1]+ '-' + 
			chem2.sequence;
	   
		outcome.prod1 = new Chemical(newseq);

		// However, we have to be consistent with the energy delta of this reaction
		var deltaE = (outcome.prod1.energy) -
					(chem1.energy + chem2.energy);
	   
		if (Math.random()<Math.exp(-deltaE/TEMPERATURE))
			return outcome;
		else return null;
	}
	
	return null;
}

function tryDblBondReact(chem1, chem2)
{
   var outcome = {
      prod1: null,
      prod2: null
   };
   
   var site = null;
   var bestdist = 0
   
   // -a=b- + -c-d-  -> -a-d- + -c-b-

   // This check is asymmetric between chem1 and chem2 because the selection of chemicals from the bath will randomly 
   for (var i=1;i<chem1.sequence.length-1;i++)
   {
		if (chem1.sequence[i]=='=') // viable site
		{
			for (var j=1;j<chem2.sequence.length-1;j++)
			{
				if (chem2.sequence[j]=='-')
				{
					var dchem=vectAdd(chem1.chemEnv[i],scalMult(-1,chem2.chemEnv[j]));
					var dist = dot(dchem,dchem);
					
					// Only choose the one best reaction as far as chemical environment
					if (dist>bestdist) 
					{
						bestdist=dist;
						site = { si: i, sj: j };
					}
				}
			}
		}
   }
   
   // I've found a place to react, so lets try to react
   if (site)
   {	   	   
		var subseq1b = chem1.sequence.slice(site.si+1, chem1.sequence.length);
		var subseq1a = chem1.sequence.slice(0,site.si);

		var subseq2a = chem2.sequence.slice(site.sj+1, chem2.sequence.length);
		var subseq2b = chem2.sequence.slice(0,site.sj);
	   
		outcome.prod1 = new Chemical(subseq1a + '-' + subseq2b);
		outcome.prod2 = new Chemical(subseq2a + '-' + subseq1b);

		// However, we have to be consistent with the energy delta of this reaction
		var deltaE = (outcome.prod1.energy + outcome.prod2.energy) -
					(chem1.energy + chem2.energy);
	   
		if (Math.random()<Math.exp(-deltaE/TEMPERATURE))
			return outcome;
		else return null;
   }
      
   // No reaction
   return null;
}

function delChem(id)
{
	chemBath[id]=chemBath[chemBath.length-1];
	chemBath.pop();
}

function doRandomReaction()
{
	var c1_id=Math.floor(Math.random()*chemBath.length);
	var c2_id;
	
	do
	{
		c2_id=Math.floor(Math.random()*chemBath.length);
	} while (c2_id == c1_id);
	
	var rtype=Math.floor(Math.random()*3);
	
	var outcome;
	
	switch (rtype)
	{
		case 0: // Decomposition
			outcome=tryDecompose(chemBath[c1_id]);
			if (outcome)
			{
				if (outcome.prod1.sequence.length>maxLength) return;
				if (outcome.prod2.sequence.length>maxLength) return;
				rxnList.push({E11: chemBath[c1_id].energy, 
					E12: -1000,
					E21: outcome.prod1.energy,
					E22: outcome.prod2.energy});
/*				console.log(chemBath[c1_id].sequence + " -> " + 
							outcome.prod1.sequence + " + " + outcome.prod2.sequence);*/
				
				chemBath[c1_id]=outcome.prod1;
				chemBath.push(outcome.prod2);
			}
			break;
			
		case 1: // Condensation
			outcome=tryRecompose(chemBath[c1_id],chemBath[c2_id]);
			if (outcome)
			{
				if (outcome.prod1.sequence.length>maxLength) return;
				rxnList.push({E11: chemBath[c1_id].energy, 
					E12: chemBath[c2_id].energy,
					E21: outcome.prod1.energy,
					E22: -1000});
/*				console.log(chemBath[c1_id].sequence + " + " + 
							chemBath[c2_id].sequence + " -> " + outcome.prod1.sequence);*/
				chemBath[c1_id]=outcome.prod1;
				
				delChem(c2_id);
			}
			break;
			
		case 2: // Exchange
			outcome=tryDblBondReact(chemBath[c1_id],chemBath[c2_id]);
			if (outcome)
			{
				if (outcome.prod1.sequence.length>maxLength) return;
				if (outcome.prod2.sequence.length>maxLength) return;
				
				rxnList.push({E11: chemBath[c1_id].energy, E12: chemBath[c2_id].energy, E21: outcome.prod1.energy, E22: outcome.prod2.energy});
/*				console.log(chemBath[c1_id].sequence + " + " + 
							chemBath[c2_id].sequence + " -> "+ outcome.prod1.sequence + " + " + outcome.prod2.sequence);*/
				
				chemBath[c1_id]=outcome.prod1;
				chemBath[c2_id]=outcome.prod2;
			}
			break;
	}
}

function initializeBath()
{
	for (var i=0;i<1500;i++)
		chemBath.push(new Chemical("B=C"));//A-B-D-D-C-B")); // "A=A-B-C=D"
}

function initializeChemistry()
{
	setupBondMatrix();
	setupBondEnergies();
	setupChemID();
	initializeBath();
}

function hashBath()
{
	var chemHash = [];
	
	for (var i=0;i<chemBath.length;i++)
	{
		if (!chemHash[chemBath[i].sequence])
		{
			chemHash[chemBath[i].sequence] = { count: 1, energy: chemBath[i].energy };
		}
		else chemHash[chemBath[i].sequence].count++;
	}
	
	return chemHash;
}

function checkMolString(molString)
{
	if (molString.length<3) return 0;
	if (molString.length%2!=1) return 0;
	
	for (var i=0;i<molString.length;i+=2)
	{
		if ((molString[i]!='A')&&(molString[i]!='B')&&(molString[i]!='C')&&(molString[i]!='D')) return 0;
	}
	for (var i=1;i<molString.length;i+=2)
	{
		if ((molString[i]!='-')&&(molString[i]!='=')) return 0;
	}
		
	return 1;
}

function addMolecule()
{
	var molString = document.getElementById("molstring").value;
	if (!checkMolString(molString)) return;
	var ammt=document.getElementById("ammt").value;
	
	for (var i=0;i<ammt;i++)
		chemBath.push(new Chemical(molString));
}

var isResetting = 0;

function resetBath()
{
	if (isResetting) return;
	
	var molString = document.getElementById("bathstring").value;
	if (!checkMolString(molString)) return;
	
	isResetting = 1;
	chemBath = [];
	
	for (var i=0;i<1500;i++)
	{
		chemBath.push(new Chemical(molString));
	}
	
	isResetting = 2;
}

function adjTemperature()
{
	TEMPERATURE = parseFloat(document.getElementById("temp").value);
}
