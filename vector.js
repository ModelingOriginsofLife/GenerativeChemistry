function normalize(v)
{
	var vn = new Vector(v.N);
	var mag = 0;
	
	for (var i=0;i<v.N;i++)
		mag += v.x[i]*v.x[i];
		
	mag = Math.sqrt(mag);
	
	for (var i=0;i<v.N;i++)
		vn.x[i] = v.x[i]/mag;
		
	return vn;
}

function Matrix(N)
{
   this.N=N;
   this.x=[];
   
   for (var j=0;j<N;j++)
      for (var i=0;i<N;i++)
         this.x[i+j*N]=(i==j);   
}

Matrix.prototype.identity = function()
{
   var N=this.N;
   
   for (var j=0;j<N;j++)
      for (var i=0;i<N;i++)
         this.x[i+j*N]=(i==j);   
}

function mulMatrix(A, B)
{
   if (A.N != B.N) return null;
   
   var C = new Matrix(A.N);
   
   for (var j=0;j<A.N;j++)
      for (var i=0;i<A.N;i++)
      {
         var val=0;
	 
	 for (var k=0;k<A.N;k++)
	 {
	    val+=A.x[k+i*A.N]*B.x[j+k*A.N];
	 }
	 
	 C.x[j+i*A.N]=val;
      }
      
   return C;
}

// matrix * vector
function matrixMulVector(A,B)
{
   if (A.N != B.N) return null;
   
   var C = new Vector(A.N);
   
   for (var i=0;i<A.N;i++)
   {
      var val=0;
	 
      for (var k=0;k<A.N;k++)
      {
         val+=A.x[k+i*A.N]*B.x[i];
      }
	 
      C.x[i]=val;
   }
   
   return C;
}

function Vector(N)
{
   this.N=N;
   this.x=[];
   for (var i=0;i<N;i++) this.x[i]=0;
}

function vectAdd(a,b)
{
   if (a.N != b.N) return null;
   var c=new Vector(a.N);
   
   for (var i=0;i<a.N;i++)
     c.x[i]=a.x[i]+b.x[i];
     
   return c;
}

function scalMult(a,b)
{
   var c=new Vector(b.N);
   
   for (var i=0;i<b.N;i++)
     c.x[i]=a*b.x[i];
     
   return c;
}

function dot(a,b)
{
   if (a.N!=b.N) return null;
   
   var c=0;
   
   for (var i=0;i<a.N;i++)
      c+=a.x[i]*b.x[i];
      
   return c;
}
