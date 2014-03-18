window.requestAnimFrame = (function(callback)
{
     return window.requestAnimationFrame || window.webkitRequestAnimationFrame || window.mozRequestAnimationFrame || window.oRequestAnimationFrame || window.msRequestAnimationFrame ||
       function(callback) {
            window.setTimeout(callback, 1000 / 30);
        };
})();

var display_element = null;
var counter=0;

function Display()
{
	var display_html = "<p>Reactions: "+counter+", ";
		
	var hash = hashBath();
	var sortHash = [];
	var total=0;
	
	for (var key in hash)
	{
		sortHash.push({key: key, value: hash[key]});
		total+=hash[key];
	}
	
	sortHash.sort(function(a,b) { return b.value - a.value; });
	
	display_html = display_html + "Total compounds: " + total + ", ";
	display_html = display_html + "Distinct compounds: " + sortHash.length + "</p>";
	
	display_html = display_html + "<ul>";
	for (var i=0;i<sortHash.length;i++)
	{
		display_html = display_html + "<li>" + sortHash[i].key + ": " + sortHash[i].value + "</li>";
	}
	
	display_html = display_html + "</ul>";
	
	if (display_element)
		display_element.innerHTML = display_html;
}

function Iterate()
{
	for (var i=0;i<100;i++)
		doRandomReaction();

	counter+=100;
	
	chemBath.push(new Chemical("D=C"));
	
	Display();
	requestAnimFrame( Iterate );
}

window.onload = function()
{
	display_element = document.getElementById("compounds");
	initializeChemistry();
	
	requestAnimFrame( Iterate );
}
