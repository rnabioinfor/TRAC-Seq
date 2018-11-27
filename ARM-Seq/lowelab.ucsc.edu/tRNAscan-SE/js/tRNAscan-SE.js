// tRNAscan-SE.js

// onchange organism
function select_org(value)
{
	if ((value == "Mitomammal") || (value == "Mitovert"))
	{
		document.getElementById("pseudogene").disabled = true;
		document.getElementById("origin").disabled = true;
//		document.getElementById("fpos").disabled = true;
		document.getElementById("breakdown").disabled = true;
		document.getElementById("gcode").value = "Vert";
	}
	else
	{
		document.getElementById("pseudogene").disabled = false;
		document.getElementById("origin").disabled = false;
//		document.getElementById("fpos").disabled = false;
		document.getElementById("breakdown").disabled = false;
		document.getElementById("gcode").value = "Universal";
	}
}