window.requestAnimFrame = (function(callback)
{
     return window.requestAnimationFrame || window.webkitRequestAnimationFrame || window.mozRequestAnimationFrame || window.oRequestAnimationFrame || window.msRequestAnimationFrame ||
       function(callback) {
            window.setTimeout(callback, 1000 / 30);
        };
})();

var compounds_element = null;
var status_element = null;
var reaction_graph = null;
var counter=0;
var context = null;

function arrowLine(x0,y0,x1,y1,size)
{
	var dx=(x1-x0), dy=(y1-y0);
	var r=Math.sqrt(dx*dx+dy*dy);
	if (r<1e-5) r=1e-5;
	dx/=r; dy/=r;
	
	var du=-dy, dv=dx;
	
	context.beginPath();
	context.moveTo(x0,y0);
	context.lineTo(x1,y1);
	context.lineTo(x1-size*dx+0.3*size*du,y1-size*dy+0.3*size*dv);
	context.moveTo(x1,y1);
	context.lineTo(x1-size*dx-0.3*size*du,y1-size*dy-0.3*size*dv);
	context.stroke();
	context.closePath();
}

function renderReactionGraph()
{
    context.fillStyle = "rgb(220,220,220)";
    context.beginPath();
    context.fillRect( 0, 0, reaction_graph.width, reaction_graph.height );
    context.closePath();
    
    context.strokeStyle = "rgb(0,0,0)";
	context.beginPath();
	context.moveTo(50,10);
	context.lineTo(50,350);
	context.lineTo(390,350);
	context.stroke();
	context.closePath();

	context.beginPath();
	context.moveTo(50,350);
	context.lineTo(390,10);
	context.stroke();
	context.closePath();
	
	context.font = "20px Arial";
    context.fillStyle = "rgb(0,0,0)";
	context.fillText("Max Energy Before", 130, 380);
	context.rotate(Math.PI/2);
	context.fillText("Max Energy After",  100, -20);
	context.rotate(-Math.PI/2);
    context.strokeStyle = "rgba(0,0,0,0.1)";
    context.fillStyle = "rgba(198,0,0,0.1)";
    
    for (var i=0;i<rxnList.length;i++)
    {
		var gx0,gy0,gx1,gy1;
		var max1,max2;
		
		max1=rxnList[i].E11; if (rxnList[i].E12>max1) max1=rxnList[i].E12;
		max2=rxnList[i].E21; if (rxnList[i].E22>max2) max2=rxnList[i].E22;
		
		gx0=50+(max1+7)*30;		
		gy0=350-(max2+7)*30;
		
		context.beginPath();
		context.arc( gx0, gy0, 3, 0, 2*Math.PI);
		context.fill();
		context.stroke();
		context.closePath();
//		console.log(gx0+","+gy0+";"+gx1+","+gy1);
		
/*		if ((gx0>=50)&&(gx0<390)&&(gy0>10)&&(gy0<=350)&&
		    (gx1>=50)&&(gx1<390)&&(gy1>10)&&(gy1<=350))
			arrowLine(gx0,gy0,gx1,gy1,10);*/
	}
}

function Display()
{
	var status_html = "<p>Reactions: "+counter+", ";
		
	var hash = hashBath();
	var sortHash = [];
	var total=0;
	
	for (var key in hash)
	{
		sortHash.push({key: key, value: hash[key].count, energy: hash[key].energy });
		total+=hash[key].count;
	}
	
	sortHash.sort(function(a,b) { return b.value - a.value; });
	
	status_html = status_html + "Total compounds: " + total + ", ";
	status_html = status_html + "Distinct compounds: " + sortHash.length + "</p>";
	
	if (status_element)
		status_element.innerHTML = status_html;
	
	var display_html = "<ul>";
	for (var i=0;i<sortHash.length;i++)
	{
		display_html = display_html + "<li>" + sortHash[i].key + ": " + sortHash[i].value + " (E = " + sortHash[i].energy.toFixed(2) + ")</li>";
	}
	
	display_html = display_html + "</ul>";
	
	if (compounds_element)
		compounds_element.innerHTML = display_html;

	renderReactionGraph();
}

function Iterate()
{
	if (rxnList.length>100)
	{
		rxnList.splice(0,rxnList.length-100);
	}
	
	for (var i=0;i<10000;i++)
		doRandomReaction();

	counter+=10000;
	
//	chemBath.push(new Chemical("D=C"));
	
	Display();
	requestAnimFrame( Iterate );
}

window.onload = function()
{
	status_element = document.getElementById("status");
	compounds_element = document.getElementById("compounds");
	reaction_graph = document.getElementById("reaction_graph");
	
	Math.seedrandom('basic_seed');
	
	initializeChemistry();
	
	context = reaction_graph.getContext('2d');
	requestAnimFrame( Iterate );
}
