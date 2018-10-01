<!DOCTYPE html>
<html>
	<head>
		<title>Processing Summary</title>
		<link 
			rel="stylesheet" 
			href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" 
			integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" 
			crossorigin="anonymous">
		<link 
			rel="stylesheet" 
			href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" 
			integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" 
			crossorigin="anonymous">
		<link 
			rel="stylesheet" 
			href="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" />
		<link 
			rel="stylesheet" 
			href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.12.1/bootstrap-table.min.css">
		<link 
			rel="stylesheet" 
			type="text/css" 
			href="https://cdn.datatables.net/1.10.18/css/dataTables.bootstrap.min.css">
		<style type="text/css">
			@-ms-viewport	 { width: device-width; }
			@-o-viewport	  { width: device-width; }
			@viewport		 { width: device-width; }
			canvas {
				-moz-user-select: none;
				-webkit-user-select: none;
				-ms-user-select: none;
			}
			.chart-container {
				width: 350px;
				margin-left: 40px;
				margin-right: 40px;
				margin-bottom: 40px;
			}
			.container {
				display: flex;
				flex-direction: row;
				flex-wrap: wrap;
				justify-content: center;
			}
		</style>
		<!--[if lt IE 9]>
			<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
			<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
		<![endif]-->
	</head>
	<body>
		<nav class="navbar navbar-inverse " role="navigation">
			<div class="container-fluid">
				<div class="navbar-header">
					<a 
						style="padding-top:4px;" 
						class="navbar-brand" 
						target="_blank" href="http://astromatic.net">
						<img 
							style="height: 50px; margin-top:0px;" 
							alt="astromatic" 
							src="http://astromatic.net/xsl/astromatic.png">
					</a>
				</div>
				<div class="nav navbar-nav navbar-right">
					<p class="navbar-text pull-right">
						<strong>
							<a target="_blank" 
								href="https://github.com/astromatic/scamp">
								<span id="soft" 
									class="label label-primary" 
									style="font-size:50;">
								</span>
							</a>
						</strong> 
						completed on <strong><span id="date"></span></strong>
						at <strong><span id="time"></span></strong> 
						using <strong><span id="nthreads"></span> </strong>
						threads (run time: <strong><span id="runtime"></span></strong>)
						started by user <strong><span id="username"></span></strong> 
						in <strong><span id="rundir"></span></strong>.
					</p>
				</div> <!-- navright -->
			</div> <!--container-fluid-->
		</nav> <!-- nav -->

		<div class="container-fluid role="main"">


			<!-- WARNINGS TABLE -->
			<div class="row-fluid">
				<div class="col-xs-12 col-sm-12 col-md-12">
					<div id ="warningDiv" class="alert alert-warning">
						<button 
							class="btn btn-warning" 
							type="button" data-toggle="collapse" 
							data-target="#warningsTableCollapse" 
							aria-expanded="false" 
							aria-controls="warningsTableCollapse">
							<strong>Warnings</strong><span id="warningsShort"></span>
						</button>
						<div id="warningsTableCollapse" class="collapse panel-body table-responsive">
							<table id="warningsTable" class="table table-responsive table-striped">
								<thead>
									<tr>
										<th>Date</th>
										<th>Time</th>
										<th>Message</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<div class="row-fluid">
				<!--  HISTOGRAMS -->
				<div class="col-xs-12 col-sm-12 col-md-12 col-lg-7">
					<div class="row">
						<div class="col-xs-12 col-sm-12 col-md-12">
							<div class="btn-group" role="group" aria-label="Basic example">
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id="plot1_btn1">Table</span>
											<span class='caret'></span>
										</button>
										<ul id="plot1_tbl" class="dropdown-menu table-dropdown">
										</ul>
									</div>
								</div>
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id ="plot1_btn2">Column</span>
											<span class='caret'></span>
										</button>
										<ul id="plot1_cols" class="dropdown-menu column-dropdown">
										</ul>
									</div>
								</div>
							</div>
							<div id="plot1">
							</div>
						</div>
					</div>
					<div class="row">
						<div class="col-xs-12 col-sm-12 col-md-12">
							<div class="btn-group" role="group" aria-label="Basic example">
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id ="plot2_btn1">Table</span>
											<span class='caret'></span>
										</button>
										<ul id="plot2_tbl" class="dropdown-menu table-dropdown">
										</ul>
									</div>
								</div>
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id ="plot2_btn2">Column</span>
											<span class='caret'></span>
										</button>
										<ul id="plot2_cols" class="dropdown-menu column-dropdown">
										</ul>
									</div>
								</div>
							</div>
							<div id="plot2">
							</div>
						</div>
					</div>
				</div> <!-- end col -->
				<!-- ALADIN -->
				<div class="col-xs-12 col-sm-12 col-md-12 col-lg-5">
					<div class="panel panel-default">
						<div class="panel-heading text-center">Catalogs Footprint</div>
						<div id="aladin-lite-div" class="panel-body" style="width: 100%; height:600px;"></div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<!-- FIELD TABLE -->
			<div class="row-fluid">
				<div class="col-md-12">
					<div id="inputFiles" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#fieldsTableCollapse" 
								aria-expanded="false" 
								aria-controls="fieldsTableCollapse">
								<strong>Summary Table on Input Files</strong> (Fields)
							</button>
						</div>
						<div id="fieldsTableCollapse" class="collapse panel-body table-responsive">
							<table 
								id="fieldsTable" 
								class="table table-hover table-bordered table-striped table-sm" 
								cellspacing="0">
								<thead>
									<tr>
										<th data-sortable="true">#</th>
										<th data-sortable="true" data-field="Filename">Filename</th>
										<th>Identifier</th>
										<th>Next</th>
										<th>Ndet</th>
										<th>Flags</th>
										<th>G</th>
										<th>A</th>
										<th>P</th>
										<th>Date</th>
										<th>Exposure Time</th>
										<th>Air Mass</th>
										<th>Right Ascension</th>
										<th>Declination</th>
										<th>Radius</th>
										<th>Pixel scale</th>
										<th class="showmatch">&#916;Pixel scale</th>
										<th class="showmatch">&#916;Position Angle</th>
										<th class="showmatch">A/S contrast</th>
										<th class="showmatch">&#916;X</th>
										<th class="showmatch">&#916;Y</th>
										<th class="showmatch">X/Y contrast</th>
										<th>&#967;2int</th>
										<th>&#967;2int High S/N</th>
										<th>&#967;2ref</th>
										<th>&#967;2ref High S/N</th>
										<th>Mag &#916;SP</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<!-- GROUPS TABLE -->
			<div class="row-fluid">
				<div class="col-md-12">
					<div id="groupProperties" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#groupPropertiesCollapse" 
								aria-expanded="false" 
								aria-controls="groupPropertiesCollapse">
								<strong>Group Properties</strong> (Fgroups)
							</button>
						</div>
					<div id="groupPropertiesCollapse" class="collapse panel-body table-responsive">
						<table id="groupsTable" class="table table-bordered table-striped">
							<thead>
								<tr>
									<th>Group name</th>
									<th class="showplot">Group Plot</th>
									<th>Index</th>
									<th>NFields</th>
									<th>Right Ascension</th>
									<th>Declination</th>
									<th>Pixel scale</th>
									<th>Maximum radius</th>
									<th>Astrom. Ref. Catalog</th>
									<th>Astrom. Ref. Band</th>
									<th class="showplot">&#967;2 Plot</th>
									<th>Astrom. &#963;int</th>
									<th>Astrom. &#961;int</th>
									<th>Astrom. &#967;2int</th>
									<th>Astrom. N2int</th>
									<th>Astrom. &#963;int High S/N</th>
									<th>Astrom. &#961;int High S/N</th>
									<th>Astrom. &#947;2int High S/N</th>
									<th>Astrom. N2int High S/N</th>
									<th class="showplot">Astrom. 1-D Int. Error Plot</th>
									<th class="showplot">Astrom. 2-D Int. Error Plot</th>
									<th>Astrom. &#916;RA ref, &#916; DEC ref</th>
									<th>Astrom. &#963;ref</th>
									<th>Astrom. &#961;ref</th>
									<th>Astrom. &#967;2ref</th>
									<th>Astrom. Nref</th>
									<th>Astrom. &#916;RA ref, &#916; DEC ref High S/N</th>
									<th>Astrom. &#963;ref High S/N</th>
									<th>Astrom. &#961;ref High S/N</th>
									<th>Astrom. &#947;2ref High S/N</th>
									<th>Astrom. Nref High S/N</th>
									<th class="showplot">Astrom. 1-D Ref. Error Plot</th>
									<th class="showplot">Astrom. 2-D Ref. Error Plot</th>
									<th>Photom. instruments</th>
									<th>Photom. &#963;int</th>
									<th>Photom. &#967;2int</th>
									<th>Photom. Nint</th>
									<th>Photom. &#963; High S/N</th>
									<th>Astrom. &#947;2int High S/N</th>
									<th>Photom. Nint High S/N</th>
									<th>Photom. &#963;ref</th>
									<th>Photom. &#967;2ref</th>
									<th>Photom. Nref</th>
									<th>Photom. &#963;ref High S/N</th>
									<th>Photom. &#947;2ref High S/N</th>
									<th>Photom. Nref High S/N</th>
									<th class="showplot">Photom. Internal Error Plot</th>
								</tr>
							</thead>
							<tbody>
							</tbody>
						</table>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<div class="row-fluid">
				<!-- ASTROINSTRUMENTS TABLE -->
				<div class="col-md-6"> 
					<div id="astroInstruments" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#astrometricInstrumentsCollapse" 
								aria-expanded="false" 
								aria-controls="astrometricInstrumentsCollapse">
								<strong>Astrometric Instruments</strong> (AstroInstruments)
							</button>
						</div>
						<div id="astrometricInstrumentsCollapse" class="collapse panel-body table-responsive">
							<table id="astrometricInstrumentsTable" class="table table-striped table-bordered" >
								<thead>
									<tr>
										<th>Name</th>
										<th>Index</th>
										<th>NFields</th>
										<th>Number of Keywords</th>
										<th>Keywords</th>
										<th class="showplot">Distortion Plot</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->

				<!-- PHOTOINSTRUMENTS TABLE -->
				<div class="col-md-6"> 
					<div id="photoInstruments" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#photoInstrumentsCollapse" 
								aria-expanded="false" 
								aria-controls="photoInstrumentsCollapse">
								<strong>Photometric Instruments</strong> (PhotInstruments)
							</button>
						</div>
						<div id="photoInstrumentsCollapse" class="collapse panel-body table-responsive">
							<table id="photometricInstrumentsTable" class="table table-striped table-bordered">
								<thead>
									<tr>
										<th>Name</th>
										<th>Index</th>
										<th>NFields</th>
										<th>Output ZP</th>
										<th>Number of Keywords</th>
										<th>Keywords</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<!-- SCAMP CONFIG -->
			<div class="row-fluid">
				<div class="col-md-12"> 
					<div id="configFile" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#configCollapse" 
								aria-expanded="false" 
								aria-controls="configCollapse">
								<strong>Scamp Configuration</strong>
							</button> 
							<span> 
								<a  href="javascript:void(0);" onclick="downloadScampConf();"> scamp.conf </a> 
							</span>
						</div>
						<div id="configCollapse" class="collapse panel-body table-responsive">
							<table id="configTable" class="table t table-responsive table-striped table-bordered table-sm">
								<thead>
									<tr>
										<th>Config Parameter</th>
										<th>Value</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->
		</div> <!-- end main container -->



		<!-- MODAL DIALOG image viewer -->
		<div 
			class="modal fade" 
			id="imageModal" 
			tabindex="-1" 
			role="dialog" 
			aria-labelledby="imageModalLabel" 
			aria-hidden="true">
			<div class="modal-dialog modal-dialog-centered" role="document" style="width:1000px;">
				<div class="modal-content">
					<div class="modal-header">
						<h5 class="modal-title" id="imageModalLabel">New message</h5>
						<button type="button" class="close" data-dismiss="modal" aria-label="Close">
							<span aria-hidden="true">&times;</span>
						</button>
					</div>
					<div class="modal-body">
						<img class="modal-img" style="margin-left:auto;margin-right:auto;display:block;"></img>
					</div>
					<div class="modal-footer">
						<button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
					</div>
				</div>
			</div>
		</div>




		<!-- SCRIPTS INCLUDES -->
		<script
			src="https://code.jquery.com/jquery-3.3.1.min.js"
			integrity="sha256-FgpCb/KJQlLNfOu91ta32o/NMZxltwRo8QtmkMRdAu8="
			crossorigin="anonymous"></script>
		<script 
			src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" 
			integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" 
			crossorigin="anonymous"></script>
		<script 
			type="text/javascript" 
			charset="utf8" 
			src="https://cdn.datatables.net/1.10.18/js/jquery.dataTables.min.js"></script>
		<script 
			type="text/javascript" 
			charset="utf8" 
			src="https://cdn.datatables.net/1.10.18/js/dataTables.bootstrap.min.js"></script>
		<script 
			type="text/javascript"
			charset="utf-8"
			src="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js"></script>
		<script 
			type="text/javascript"
			charset="utf-8"
			src="https://cdn.plot.ly/plotly-1.41.0.min.js"></script>
		<script>
			/*! 
			 * IE10 viewport hack for Surface/desktop Windows 8 bug
			 * Copyright 2014-2015 Twitter, Inc.
			 * Licensed under MIT (https://github.com/twbs/bootstrap/blob/master/LICENSE)
			 */
			(
				function () {
					'use strict';
					if (navigator.userAgent.match(/IEMobile\/10\.0/)) {
						var msViewportStyle = document.createElement('style')
						msViewportStyle.appendChild(
							document.createTextNode(
								'@-ms-viewport{width:auto!important}'
							)
						)
						document.querySelector('head').appendChild(msViewportStyle)
					}
				} 
			)();
		</script>
	</body>
</html>

<script>
/*
 * Data contained in "scamp_data" object, is loaded at the very end 
 * of this html file is his own <script> tag.
 */

/* from an array of object (data), return the "value" property of   
 * object havinig "str" as "name" property. */   
function getElemVal(str, data) {
	var value = "";
	$.each(data, function (i, elem) {
			if (elem.name == str) {
			value = elem.value;
			return;
			}
			});
	return value;
}

function getFlagValHelper(b, t) {
	if (b)
		return t;
	return "-";
}

function getElemAverageValHelper(a) {
	var value = 0.0;
	var n = 0;
	$.each(a, function(i, elem) {
			value += elem;
			n++;
			});
	return value / n;
}

function getElemListValHelperFixedHelper(arr, unit, fix) {
	var value = "";
	$.each(arr, function(i, elem) {
			if (fix >= 0) {
			elem = elem.toFixed(fix);
			}
			value += elem + unit + " ";
			});
	return value;
}

function getElemListValHelper(arr, unit) {
	return getElemListValHelperFixedHelper(arr, unit, -1);
}

function getRaValHelper(value) {
	var a = Math.floor(value[0] / 15.0);
	var b = Math.floor((value[0] * 4) % 60);
	var c = Math.floor((value[0] *240) % 60);
	return a + ":" + b + ":" + c.toFixed(2);
}

function getDecValHelper(value) {
	var sign = "";
	if (value[1] < 0) {
		sign = "-";
		value[1] = 0 - value[1];
	}
	var a = Math.floor(value[1]);
	var b = Math.floor((value[1] * 60) % 60);
	var c = Math.floor((value[1] * 3600) % 60);

	return sign + a + ":" + b + ":" + c.toFixed(2);
}

function setAladinPos() {
	var left_max	= 1000;
	var top_max	 = -1;
	var right_max   = -1;
	var bottom_max  = 1000;
	$.each(scamp_data.Fields, function(i, field) {
		for (var i=0; i<field.Set_Polygon.value.length; i++) {
			var poly = field.Set_Polygon.value[i];
			if (poly[0][0] < left_max)
				left_max = poly[0][0];
			if (poly[0][1] > top_max)
				top_max = poly[0][1];
			if (poly[1][0] > right_max)
				right_max = poly[1][0];
			if (poly[2][1] < bottom_max)
				bottom_max = poly[2][1];
		}
	})

	var center_lng = (left_max + right_max) / 2;
	var center_lat = (top_max + bottom_max) / 2;
	aladin.gotoRaDec(center_lng, center_lat);

	var height_max = Math.abs(top_max - bottom_max);
	aladin.setFov(height_max * 4);
}


PLOT_COLORS = {
	'#plot1': '#fce84f',
	'#plot2': '#4e9a06'
};

FIRST_COL = {
	'Fields': 'Observation_Date',
	'Fgroups': 'Max_Radius',
	'AstroInstruments': 'NFields',
	'PhotInstruments': 'NFields'
}

function setDDmenuCallback() {
	$(".dropdown-menu li a").click(function() {
		var value   = $(this).text();
		var listId  = $(this).parent().parent().attr('id');
		var isColumn = listId.endsWith('_cols');
		var divId = '';
		if (isColumn) {
			divId = listId.substring(0, listId.length - 5)
			var tableId = $('#' + divId + '_btn1').text();
			generateHistogram(tableId, value, '#' + divId, PLOT_COLORS['#'+divId]);
		} else {
			divId = listId.substring(0, listId.length - 4)
			generateHistogram(value, FIRST_COL[value], '#' + divId, PLOT_COLORS['#'+divId]);
		}
	});
}


COL_Fields = `
<li class="dropdown-header">Column</li>
<li role="separator" class="divider"></li>
<li><a href="javascript:void(0);">Observation_Date</a></li>
<li><a href="javascript:void(0);">Exposure_Time</a></li>
<li><a href="javascript:void(0);">Air_Mass</a></li>
<li><a href="javascript:void(0);">Pixel_Scale</a></li>
<li><a href="javascript:void(0);">AS_Contrast</a></li>
<li><a href="javascript:void(0);">XY_Contrast</a></li>
`;

COL_Groups = `
<li class="dropdown-header">Column</li>
<li role="separator" class="divider"></li>
<li><a href="javascript:void(0);">Max_Radius</a></li>
<li><a href="javascript:void(0);">AstromNDets_Internal</a></li>
`;

COL_AstroInstruments = `
<li class="dropdown-header">Column</li>
<li role="separator" class="divider"></li>
<li><a href="javascript:void(0);">NFields</a></li>
`

COL_PhotInstruments = `
<li class="dropdown-header">Column</li>
<li role="separator" class="divider"></li>
<li><a href="javascript:void(0);">NFields</a></li>
`

function generateHistogram(table, column, divId, color) {
	$(divId + '_btn1').html(table);
	$(divId + '_btn2').html(column);

	if (table == 'Fields') {
		$(divId + '_cols').html(COL_Fields);
	} else if (table == 'Fgroups') {
		$(divId + '_cols').html(COL_Groups);
	} else if (table == 'PhotInstruments') {
		$(divId + '_cols').html(COL_PhotInstruments);
	} else if (table == 'AstroInstruments') {
		$(divId + '_cols').html(COL_AstroInstruments);
	}
	setDDmenuCallback();

	var d3 = Plotly.d3;
	var gd = d3.select(divId).node();
	var x = [];
	$.each(scamp_data[table], function(i, row) {
		x[i] = row[column].value;
	});

	var layout = {
		showlegend:  false,
		xaxis: {title: column},
		margin: {l:50, r:50, b:80, t:40, pad:4},
		height: 300
	};

	Plotly.newPlot(gd, [{
		type: 'histogram',
		autobinx: false,
		histfunc: "density",
		x: x,
		marker: {color: color, line: {color: '#555',width:1}},
		opacity: 0.8,
	}], layout, {displayModeBar:true,showLink:false,displayLogo:false});

	return gd;
}

function aladinDraw(higlight) {
	aladin.removeLayers();
	var lineWidth = 1;
	var color = '#bbbbbbff';
	var overlay;
	var poly;
	var drawline;

	$.each(scamp_data.Fields, function(i, field) {
		if (higlight == i)
			return;

		overlay = A.graphicOverlay({color: color, lineWidth: lineWidth});
		aladin.addOverlay(overlay);

		for (var i=0; i<field.Set_Polygon.value.length; i++) {
			poly = field.Set_Polygon.value[i];
			drawline = poly;
			drawline.push(poly[0]);
			overlay.add(A.polyline(drawline));
		}
	})

	if (higlight < 0)
		return;

	/* remain the hiligted one */
	lineWidth = 10;
	color = '#ff0000ff';
	var field = scamp_data.Fields[higlight];
	var overlay = A.graphicOverlay({color: color, lineWidth: lineWidth});
	aladin.addOverlay(overlay);
	for (var i=0; i<field.Set_Polygon.value.length; i++) {
		var poly = field.Set_Polygon.value[i];
		var drawline = poly;
		drawline.push(poly[0]);
		overlay.add(A.polyline(drawline));
	}
}

function generateImageColHelper(imageUrl) {
	if (imageUrl == null)
		return "<td></td>";

	var value = "";
	/* value += "<td><a type='button' rel='popover' data-img='"+imageUrl+"'>"; */
	value += "<td><a type='button' data-toggle='modal' data-target='#imageModal' data-imgurl='" + imageUrl.value + "'>";
	value += "<img width='100' class='img-fluid' src='" + imageUrl.value + "' />";
	value += "</a></td>";
	return value;
}

function downloadScampConf() {
	/* generate scamp conf string */
	var text = "";
	for (var i = 0; i < scamp_data.Configuration.length; i++) {
		var configObj = scamp_data.Configuration[i];
		text += configObj.name + "\t" + configObj.value + '\n';
	}
	var element = document.createElement('a');
	element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
	element.setAttribute('download', 'scamp.conf');
	element.style.display = 'none';
	document.body.appendChild(element);
	element.click();
	document.body.removeChild(element);
}

$(document).ready(function() {
	/* first initialize aladin */
	aladin = A.aladin('#aladin-lite-div', {
		survey: "P/DSS2/color", 
		fov: 60,
		showReticle: true,
		showZoomControl: true,
		showFullScreenControl: true,
		showLayersControl: true,
		showGotoControl: true,
		showShareControl: true,
		fulscreen: true
	});

	aladinDraw(-1);
	console.log(scamp_data);

	/* build status string */
	$('#soft').text(scamp_data.Software.Name.value +" "+scamp_data.Software.Version.value);
	$('#date').text(scamp_data.Software.Date.value);
	$('#time').text(scamp_data.Software.Time.value);
	$('#nthreads').text(scamp_data.Software.NThreads.value);
	$('#runtime').text(scamp_data.Software.Duration.value + " s");
	$('#username').text(scamp_data.Software.User.value);
	$('#rundir').text(scamp_data.Software.Path.value);

	/* show/hide match option and plots */
	var showmatch = getElemVal("MATCH", scamp_data.Configuration);
	var showplot  = getElemVal("CHECKPLOT_DEV", scamp_data.Configuration)[0];
	showplot = (showplot == "PNG") ? true : false;

	/* 
	 * build fields table 
	 */
	$.each(scamp_data.Fields, function(i, field) {
		var tdClass = "<td>";
		var magLimit = 2.0;
		var fieldNum = i + 1;
		if (field.XY_Contrast.value.toFixed(1) < magLimit) {
			tdClass = '<td class="danger">';
			scamp_data.Warnings.push({
				Date: {value: scamp_data.Software.Date.value},
				Text: {value: "Catalog number " + fieldNum + " failled the crossmatch with the external catalog. (catalog file: " + field.Catalog_Name.value + ")"},
				Time: {value: scamp_data.Software.Time.value}
			});
		} 
		var table_row = "";
		table_row += '<tr>';
		table_row += tdClass +  field.Catalog_Number.value + "</td>";
		table_row += tdClass +  field.Catalog_Name.value + "</td>";
		table_row += tdClass +  field.Image_Ident.value + "</td>";
		table_row += tdClass +  field.NExtensions.value + "</td>";
		table_row += tdClass +  field.NDetect.value + "</td>";
		table_row += tdClass +  getFlagValHelper(field.Ext_Header.value, "H")  + getFlagValHelper(field.Photom_Flag, "P") +  "</td>";
		table_row += tdClass +  field.Group.value + "</td>";
		table_row += tdClass +  field.Astr_Instrum.value + "</td>";
		table_row += tdClass +  field.Phot_Instrum.value + "</td>";
		table_row += tdClass +  field.Observation_Date.value + "</td>";
		table_row += tdClass +  field.Exposure_Time.value.toFixed(3) + "</td>";
		table_row += tdClass +  field.Air_Mass.value.toFixed(2) + "</td>";
		table_row += tdClass +  getRaValHelper(field.Field_Coordinates.value) + "</td>";
		table_row += tdClass +  getDecValHelper(field.Field_Coordinates.value) + "</td>";
		table_row += tdClass +  field.Max_Radius.value.toFixed(3) + "'" + "</td>";
		table_row += tdClass +  getElemAverageValHelper(field.Pixel_Scale.value).toFixed(4) + "''" + "</td>";
		if (showmatch) {
			table_row += tdClass +  field.DPixel_Scale.value.toFixed(4) + "</td>";
			table_row += tdClass +  field.DPos_Angle.value + "°" + "</td>";
			table_row += tdClass +  field.AS_Contrast.value.toFixed(1) + "</td>";
			table_row += tdClass +  field.DX.value.toExponential() + "°" + "</td>";
			table_row += tdClass +  field.DY.value.toExponential() + "°" + "</td>";
			table_row += tdClass +  field.XY_Contrast.value.toFixed(1) + "</td>";
		} else {
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
		}
		table_row += tdClass +  field.Chi2_Internal.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.Chi2_Internal_HighSN.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.Chi2_Reference.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.Chi2_Reference_HighSN.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.ZeroPoint_Corr.value.toFixed(3) + "</td>";
		table_row += "</tr>";
		$(table_row).appendTo("#fieldsTable tbody");
	});

	$('#fieldsTable').DataTable({
		paging: false,
		info: false,
		searching: true
	});

	setAladinPos();
	$('#fieldsTable tbody tr').on('click', function(event) {
		if ($(this).hasClass('info'))  {
			$(this).removeClass('info');
			aladinDraw(-1);
			return;
		}
		$(this).addClass('info').siblings().removeClass('info');
		aladinDraw($(this).children('td:first').text() - 1);
	});

	/* 
	 * build fields groups table 
	 */
	$.each(scamp_data.Fgroups, function(i, group) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  group.Name.value + "</td>";
		if (showplot) {
			table_row += generateImageColHelper(group.FgroupsPlot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "<td>" +  group.Index.value + "</td>";
		table_row += "<td>" +  group.NFields.value + "</td>";
		table_row += "<td>" +  getRaValHelper(group.Field_Coordinates.value) + "</td>";
		table_row += "<td>" +  getDecValHelper(group.Field_Coordinates.value) + "</td>";
		table_row += "<td>" +  getElemAverageValHelper(group.Pixel_Scale.value).toFixed(4) + "''" + "</td>";
		table_row += "<td>" +  group.Max_Radius.value.toFixed(3) + "'" + "</td>";
		table_row += "<td>" +  group.AstRef_Catalog.value + "</td>";
		table_row += "<td>" +  group.AstRef_Band.value + "</td>";
		if (showplot) {
			table_row += generateImageColHelper(group.Chi2Plot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.AstromSigma_Internal.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Internal.value.toFixed(5) + "</td>";
		table_row += "<td>" +  group.AstromChi2_Internal.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Internal.value + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.AstromSigma_Internal_HighSN.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Internal_HighSN.value.toFixed(5)  + "</td>";
		table_row += "<td>" +  group.AstromChi2_Internal_HighSN.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Internal_HighSN.value + "</td>";
		if (showplot) {
			table_row += generateImageColHelper(group.IntErr1DimPlot);
			table_row += generateImageColHelper(group.IntErr2DimPlot);
		} else {
			table_row += "<td></td>";
			table_row += "<td></td>";
		}
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.AstromOffset_Reference.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.AstromSigma_Reference.value, "'' ", 3) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Reference.value.toFixed(4) + "</td>";
		table_row += "<td>" +  group.AstromChi2_Reference.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Reference.value + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.AstromOffset_Reference_HighSN.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.AstromSigma_Reference_HighSN.value, "'' ", 3) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Reference_HighSN.value.toFixed(4) + "</td>";
		table_row += "<td>" +  group.AstromChi2_Reference_HighSN.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Reference_HighSN.value + "</td>";
		if (showplot) {
			table_row += generateImageColHelper(group.RefErr1DimPlot);
			table_row += generateImageColHelper(group.RefErr2DimPlot);
		} else {
			table_row += "<td></td>";
			table_row += "<td></td>";
		}
		table_row += "<td>" +  getElemListValHelper(group.PhotInstru_Name.value,", ") + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotSigma_Internal.value, " ", 6) + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotChi2_Internal.value, " ", 4) + "</td>";
		table_row += "<td>" +  getElemListValHelper(group.PhotNDets_Internal.value, " ") + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotSigma_Internal_HighSN.value, " ", 6) + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotChi2_Internal_HighSN.value, " ", 2) + "</td>";
		table_row += "<td>" +  getElemListValHelper(group.PhotNDets_Internal_HighSN.value, " ") + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotSigma_Reference.value, " ", 6) + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotChi2_Reference.value, " ", 6) + "</td>";
		table_row += "<td>" +  getElemListValHelper(group.PhotNDets_Reference.value, " ") + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotSigma_Reference_HighSN.value, " ", 6) + "</td>";
		table_row += "<td>" +  getElemListValHelperFixedHelper(group.PhotChi2_Reference_HighSN.value, " ", 6) + "</td>";
		table_row += "<td>" +  getElemListValHelper(group.PhotNDets_Reference_HighSN.value, " ") + "</td>";
		if (showplot) {
			table_row += generateImageColHelper(group.PhotErrPlot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "</tr>";
		$(table_row).appendTo("#groupsTable tbody");
	});

	$('#groupsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});

	/* 
	 * build astrometric instruments table 
	 */
	$.each(scamp_data.AstroInstruments, function(i, astroinstru) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  astroinstru.Name.value + "</td>";
		table_row += "<td>" +  astroinstru.Index.value + "</td>";
		table_row += "<td>" +  astroinstru.NFields.value + "</td>";
		table_row += "<td>" +  astroinstru.NKeys.value + "</td>";
		table_row += "<td>" +  getElemListValHelper(astroinstru.Keys.value, " ") + "</td>";
		if (showplot) {
			table_row += generateImageColHelper(astroinstru.DistPlot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "</tr>";
		$(table_row).appendTo("#astrometricInstrumentsTable tbody");
	});

	$('#astrometricInstrumentsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});


	/* 
	 * build photometric instruments table 
	 */
	$.each(scamp_data.PhotInstruments, function(i, photoinstru) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  photoinstru.Name.value + "</td>";
		table_row += "<td>" +  photoinstru.Index.value + "</td>";
		table_row += "<td>" +  photoinstru.NFields.value + "</td>";
		table_row += "<td>" +  photoinstru.MagZeroPoint_Output.value + "</td>";
		table_row += "<td>" +  photoinstru.NKeys.value + "</td>";
		table_row += "<td>" +  getElemListValHelper(photoinstru.Keys.value, " ") + "</td>";
		table_row += "</tr>";
		$(table_row).appendTo("#photometricInstrumentsTable tbody");
	});

	$('#photometricInstrumentsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});


	/* 
	 * build configuration table 
	 */
	$.each(scamp_data.Configuration, function(i, config) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  config.name + "</td>";

		var value = "";
		if (config.datatype.includes("array")) {
			for (var i = 0; i < config.value.length; i++) {
				value += config.value[i] + ", ";
			}
		} else {
			value = config.value;
		}
		table_row += "<td>" +  value + "</td>";
		table_row += "</tr>";
		$(table_row).appendTo("#configTable tbody");
	});

	$('#configTable').DataTable({
		scrollY: "400px",
		scrollCollapse: true,
		paging: false,
		info: false,
		searching: false
	});

	/* 
	 * build warnings table 
	 */
	$.each(scamp_data.Warnings, function(i, warn) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  warn.Date.value + "</td>";
		table_row += "<td>" +  warn.Time.value + "</td>";
		table_row += "<td>" +  warn.Text.value + "</td>";
		$(table_row).appendTo("#warningsTable tbody");
	});

	$('#warningsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});

	/*
	 * Hide unused columns
	 */
	if (!showmatch) {
		$('#fieldsTable th:nth-child(17)').hide();
		$('#fieldsTable td:nth-child(17)').hide();
		$('#fieldsTable th:nth-child(18)').hide();
		$('#fieldsTable td:nth-child(18)').hide();
		$('#fieldsTable th:nth-child(19)').hide();
		$('#fieldsTable td:nth-child(19)').hide();
		$('#fieldsTable th:nth-child(20)').hide();
		$('#fieldsTable td:nth-child(20)').hide();
		$('#fieldsTable th:nth-child(21)').hide();
		$('#fieldsTable td:nth-child(21)').hide();
		$('#fieldsTable th:nth-child(22)').hide();
		$('#fieldsTable td:nth-child(22)').hide();
	}
	if (!showplot) {
		$('#groupsTable th:nth-child(2)').hide();
		$('#groupsTable td:nth-child(2)').hide();
		$('#groupsTable th:nth-child(11)').hide();
		$('#groupsTable td:nth-child(11)').hide();
		$('#groupsTable th:nth-child(20)').hide();
		$('#groupsTable td:nth-child(20)').hide();
		$('#groupsTable th:nth-child(21)').hide();
		$('#groupsTable td:nth-child(21)').hide();
		$('#groupsTable th:nth-child(32)').hide();
		$('#groupsTable td:nth-child(32)').hide();
		$('#groupsTable th:nth-child(33)').hide();
		$('#groupsTable td:nth-child(33)').hide();
		$('#groupsTable th:nth-child(47)').hide();
		$('#groupsTable td:nth-child(47)').hide();
		$('#astrometricInstrumentsTable th:nth-child(6)').hide();   
		$('#astrometricInstrumentsTable td:nth-child(6)').hide();  
	}

	if (scamp_data.Warnings.length == 0) {
		$("#warningDiv").hide();
	} else {
		var nwarn = $('#warningsTable tbody tr').length;
		$('#warningsShort').html(' (' + nwarn + ')');
	}


	/*
	 * Configure our modal image viewers
	 */
	$('#imageModal').on('show.bs.modal', function (event) {
		var button = $(event.relatedTarget); // Button that triggered the modal
		var img = button.data('imgurl'); // Extract info from data-* attributes
		var modal = $(this);
		modal.find('.modal-title').text("Image: " + img);
		imgtag = modal.find('.modal-img');
		imgtag.attr("src", img);
		imgtag.attr("title", img);
			imgtag.attr("alt", img);
	})

	$('a[rel=popover]').popover({
		animation: true, 
		container: 'body', 
		html: true, 
		placement: 'bottom', 
		content: function() {
			return "<img src='"+$(this).data('img') + "' />";
		}
	});


	var gd1 = generateHistogram('Fields', 'Observation_Date', '#plot1', PLOT_COLORS['#plot1']);
	var gd2 = generateHistogram('Fields', 'XY_Contrast',	  '#plot2', PLOT_COLORS['#plot2']);

	window.onresize = function() {
		Plotly.Plots.resize(gd1);
		Plotly.Plots.resize(gd2);
	}


	$(".table-dropdown").each(function(index) {
		$(this).append(`
			<li class="dropdown-header">Tables</li>
			<li role="separator" class="divider"></li>
			<li><a href="javascript:void(0);">Fields</a></li>
			<li><a href="javascript:void(0);">Fgroups</a></li>
			<li><a href="javascript:void(0);">PhotInstruments</a></li>
			<li><a href="javascript:void(0);">AstroInstruments</a></li>
		`);
	});

	setDDmenuCallback();

	$("#configCollapse").on('shown.bs.collapse', function() {
		$('#configTable').DataTable().columns.adjust().draw();
	});

});
</script>
