<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
        <!ENTITY Delta "&#916;">
        <!ENTITY chi "&#967;">
        <!ENTITY br "&#60;&#98;&#115;&#32;&#47;&#62;">
	]>

<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- *********************** Global XSL template ************************** -->
 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <html>

<!-- HTML head  -->

   <head>

<!-- javascript -->

<!--  <script type="text/javascript" language="javascript"> -->

    <script src="http://terapix.iap.fr/cplt/xsl/sorttable.js"/>

    <style type="text/css">
     p.sansserif {font-family: sans-serif}
     p.italic {font-style: italic}
     body {background-color: white}
     mono {font-family: monospace}
     elen {font-family: monospace; font-size: 100%; font-weight: bold; color: green }
     elep {font-family: monospace; font-size: 100%; font-weight: bold; color: red }
     el {font-family: monospace; font-size: 100%; color: black}
     a {text-decoration: none}

     table.sortable a.sortheader
      {
      background-color:#FFEECC;
      color: black;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }
     table.sortable span.sortarrow
      {
      color: black;
      font-weight: bold;
      text-decoration: none;
      }
     table.sortable a.sortheader.sub
      {
      vertical-align: sub;
      }
    </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </head>

<!-- HTML body -->

    <BODY>
     <TABLE BORDER="0" CELLPADDING="0" CELLSPACING="0" WIDTH="100%">
      <TR>
       <TD ALIGN="LEFT">
        <TABLE BORDER="0">
         <TR>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixLogo.png" ALT="Terapix"/>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixTitle.png" ALT="Logo"/>
          </TD>
          <TD ALIGN="CENTER" width="2000">
           <FONT color="#669933">
            <B> Processing summary</B>
           </FONT>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixPicture.gif" ALT="Terapix banner"/>
          </TD>
         </TR>
        </TABLE>
       </TD>
      </TR>
      <TR>
       <TD>
        <TABLE BORDER="0" WIDTH="100%" BGCOLOR="#000000">
         <TR>
          <TH BGCOLOR="#000000" ALIGN="LEFT"><FONT SIZE="-1" COLOR="#FFFFFF"> Home > Tools > Data reduction</FONT></TH>
         </TR>
        </TABLE>
       </TD>
      </TR>
     </TABLE>
    <xsl:call-template name="VOTable"/>
   </BODY>
  </html>
 </xsl:template>

<!-- **************** Generic XSL template for VOTables ****************** -->
 <xsl:template name="VOTable">
  <xsl:for-each select="/VOTABLE">
   <xsl:call-template name="Resource"/>
  </xsl:for-each>
 </xsl:template>

<!-- *************** Generic XSL template for Resources ****************** -->
 <xsl:template name="Resource">
  <xsl:for-each select="RESOURCE">
   <xsl:choose>
    <xsl:when test="@ID='SCAMP'">
     <xsl:call-template name="SCAMP"/>
    </xsl:when>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>

<!-- ********************** XSL template for SCAMP *********************** -->
 <xsl:template name="SCAMP">
  <xsl:for-each select="RESOURCE[@ID='MetaData']">
   <xsl:call-template name="RunInfo"/>
   <xsl:for-each select="TABLE[@ID='Fields']">
    <xsl:call-template name="Fields"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='FGroups']">
    <xsl:call-template name="FGroups"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Astrometric_Instruments']">
    <xsl:call-template name="Astr_Instr"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Photometric_Instruments']">
    <xsl:call-template name="Phot_Instr"/>
   </xsl:for-each>
   <xsl:for-each select="RESOURCE[@ID='Config']">
    <xsl:call-template name="Config"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Warnings']">
    <xsl:call-template name="Warnings"/>
   </xsl:for-each>
  </xsl:for-each>
 </xsl:template>

<!-- ************* Generic XSL RunInfo template for MetaData ************* -->
 <xsl:template name="RunInfo">
  <p>
<!-- Software name, version, date, time and number of threads -->
   <a>
    <xsl:attribute name="href">
     <xsl:value-of select="PARAM[@name='Soft_URL']/@value"/>
    </xsl:attribute>
    <b>
     <xsl:value-of select="PARAM[@name='Software']/@value"/>&nbsp;<xsl:value-of select="PARAM[@name='Version']/@value"/>
    </b>
   </a>
   started on
   <b><xsl:value-of select="PARAM[@name='Date']/@value"/></b>
   at
   <b><xsl:value-of select="PARAM[@name='Time']/@value"/></b>
   with
   <b><xsl:value-of select="PARAM[@name='NThreads']/@value"/></b>
   thread<xsl:if test="PARAM[@name='NThreads']/@value &gt; 1">s</xsl:if>

<!-- Run time -->
   <xsl:variable name="duration" select="PARAM[@name='Duration']/@value"/>
   (run time:
    <b>
     <xsl:choose> 
      <xsl:when test="$duration &gt; 3600.0">
       <xsl:value-of
	select='concat(string(floor($duration div 3600)),
	" h ", format-number(floor(($duration div 60) mod 60.0), "00"),
	" min")'/>
      </xsl:when>
      <xsl:otherwise>
       <xsl:choose>
        <xsl:when test="$duration &gt; 60.0">
         <xsl:value-of
	  select='concat(format-number(floor($duration div 60),"##"),
	  " min ", format-number(floor($duration mod 60.0), "00")," s")'/>
        </xsl:when>
        <xsl:otherwise>
         <xsl:value-of select='concat(string($duration), " s")'/>
        </xsl:otherwise>
       </xsl:choose>
      </xsl:otherwise>
     </xsl:choose>
    </b>)
    <br />
   by user <b><xsl:value-of select="PARAM[@name='User']/@value"/></b>
   from <b><xsl:value-of select="PARAM[@name='Host']/@value"/></b>
   in <b><mono><xsl:value-of select="PARAM[@name='Path']/@value"/></mono></b>
  </p>
  <p>
   <b style="color: red"><xsl:if test="PARAM[@name='Error_Msg']/@value &gt; 0">
    An Error occured!!! </xsl:if>
   <xsl:value-of select="PARAM[@name='Error_Msg']/@value"/></b>
  </p>
  <p>
  <i>click to expand or hide tables</i>
  </p>
 </xsl:template>

<!-- ********************** XSL template for Fields ********************** -->
  <xsl:template name="Fields">
   <xsl:variable name="asplotflag" select="count(PARAM[@name='AllSkyPlot'])"/>
   <xsl:variable name="name" select="count(FIELD[@name='Catalog_Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ident" select="count(FIELD[@name='Image_Ident']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="next" select="count(FIELD[@name='NExtensions']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ndet" select="count(FIELD[@name='NDetect']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="headflag" select="count(FIELD[@name='Ext_Header']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="photflag" select="count(FIELD[@name='Photom_Flag']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="group" select="count(FIELD[@name='Group']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astrinstru" select="count(FIELD[@name='Astr_Instrum']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="photinstru" select="count(FIELD[@name='Phot_Instrum']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="coord" select="count(FIELD[@name='Field_Coordinates']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="radius" select="count(FIELD[@name='Max_Radius']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pixscale" select="count(FIELD[@name='Pixel_Scale']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="dpscale" select="count(FIELD[@name='DPixelScale']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="dposangle" select="count(FIELD[@name='DPosAngle']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ascont" select="count(FIELD[@name='AS_Contrast']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="dx" select="count(FIELD[@name='DX']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="dy" select="count(FIELD[@name='DY']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="xycont" select="count(FIELD[@name='XY_Contrast']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2ref" select="count(FIELD[@name='Chi2_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="shear" select="count(FIELD[@name='Shear']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="sposangle" select="count(FIELD[@name='Shear_PosAngle']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="zpcorr" select="count(FIELD[@name='ZeroPoint_Corr']/preceding-sibling::FIELD)+1"/>
   <xsl:if test="$asplotflag &gt; 0">
    <a target="_blank">
     <xsl:attribute name="href">
      <xsl:value-of select="PARAM[@name='AllSkyPlot']/@value"/>
     </xsl:attribute>
     <img width="400">
      <xsl:attribute name="src">
       <xsl:value-of select="PARAM[@name='AllSkyPlot']/@value"/>
      </xsl:attribute>
      <xsl:attribute name="title">
       <xsl:value-of select="PARAM[@name='AllSkyPlot']/@value"/>
      </xsl:attribute>
     </img>
    </a>
   </xsl:if>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('scamp')">
     Summary Table on Input Files
    </BUTTON>
    <TABLE class="sortable" id="scamp" BORDER="2" style="display: none">
<!--
 <TR>
  <TH COLSPAN="3" BGCOLOR="#DDCC99" ALIGN="LEFT">
  <FONT SIZE="-1">
   <xsl:value-of select="ASTROM/nfield"/>
    input<xsl:if test="ASTROM/nfield &gt; 1">s</xsl:if>
  </FONT>
  </TH>
 </TR>
-->
     <TR>
      <TH BGCOLOR="#FFEECC">Filename</TH>
      <TH BGCOLOR="#FFEECC">Identifier</TH>
      <TH BGCOLOR="#FFEECC">Next</TH>
      <TH BGCOLOR="#FFEECC">Ndet</TH>
      <TH BGCOLOR="#FFEECC">Flags</TH>
      <TH BGCOLOR="#FFEECC">G</TH>
      <TH BGCOLOR="#FFEECC">A</TH>
      <TH BGCOLOR="#FFEECC">P</TH>
      <TH BGCOLOR="#FFEECC">alpha</TH>
      <TH BGCOLOR="#FFEECC">delta</TH>
      <TH BGCOLOR="#FFEECC">Radius</TH>
      <TH BGCOLOR="#FFEECC">Pixel scale</TH>
      <TH BGCOLOR="#FFEECC">&Delta;Pixel Scale</TH>
      <TH BGCOLOR="#FFEECC">&Delta;Position Angle</TH>
      <TH BGCOLOR="#FFEECC">A/S contrast</TH>
      <TH BGCOLOR="#FFEECC">&Delta;X</TH>
      <TH BGCOLOR="#FFEECC">&Delta;Y</TH>
      <TH BGCOLOR="#FFEECC">X/Y contrast</TH>
      <TH BGCOLOR="#FFEECC">&chi;<sup>2</sup> <i>w.r.t.</i> Reference Astrometric Catalog</TH>
      <TH BGCOLOR="#FFEECC">Shear</TH>
      <TH BGCOLOR="#FFEECC">Shear Position Angle</TH>
      <TH BGCOLOR="#FFEECC">MagZP.corr</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ident]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$next]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ndet]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$headflag] = 'T'">
           <elen>H</elen>
          </xsl:when>
          <xsl:otherwise>
            <el>-</el>
          </xsl:otherwise>
         </xsl:choose>
         <xsl:choose>
          <xsl:when test="TD[$photflag] = 'T'">
            <elen>P</elen>
          </xsl:when>
          <xsl:otherwise>
            <el>-</el>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$group]"/></el>
        </td>
        <td align="left" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astrinstru]"/></el>
        </td>
        <td align="left" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$photinstru]"/></el>
        </td>
<!-- Alpha -->
        <td align="center" BGCOLOR="#EEEEEE">
         <el>
          <xsl:variable name="alpha" select="substring-before(string(TD[$coord]), ' ')"/>
          <xsl:value-of
		select='concat(format-number(floor($alpha div 15.0), "00"),":",
		format-number(floor(($alpha * 4) mod 60.0), "00"),":",
		format-number(floor(($alpha * 240.0) mod 60.0), "00.00"))'/>
         </el>
        </td>
<!-- Delta -->
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:variable name="delta" select="substring-after(string(TD[$coord]), ' ')"/>
         <el>
          <xsl:choose>
           <xsl:when test="$delta &lt; 0.0">
            <xsl:value-of
		select='concat("-", format-number(floor(-$delta), "00"),":",
		format-number(floor((-$delta * 60) mod 60.0), "00"),":",
		format-number(floor((-$delta * 3600.0) mod 60.0), "00.0"))'/>
           </xsl:when>
           <xsl:otherwise>
            <xsl:value-of
		select='concat("+", format-number(floor($delta), "00"),":",
		format-number(floor(($delta * 60) mod 60.0), "00"),":",
		format-number(floor(($delta * 3600.0) mod 60.0), "00.0"))'/>
           </xsl:otherwise>
          </xsl:choose>
         </el>
        </td>
<!-- Radius -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$radius],'##0.000')"/>&amin;</el>
        </td>
<!-- Pixel scale --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="pix1" select="number(substring-before(string(TD[$pixscale]), ' '))"/>
         <xsl:variable name="pix2" select="number(substring-after(string(TD[$pixscale]), ' '))"/>
         <el><xsl:value-of select="format-number(($pix1+$pix2) div 2.0, '##0.0000')"/>&asec;</el>
        </td>
<!-- Delta Pixel scale --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$dpscale], '##0.0000')"/></el>
        </td>
<!-- Delta Pos Angle --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$dposangle]"/>&deg;</el>
        </td>
<!-- A/S contrast --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$xycont] &lt; 2.0">
           <elep><xsl:value-of select="format-number(TD[$ascont], '##0.0')"/></elep>
          </xsl:when>
          <xsl:otherwise>
            <elen><xsl:value-of select="format-number(TD[$ascont], '##0.0')"/></elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
<!-- Delta X -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$dx]"/>&deg;</el>
        </td>
<!-- Delta Y -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$dy]"/>&deg;</el>
        </td>
<!-- X/Y contrast --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$xycont] &lt; 2.0">
           <elep><xsl:value-of select="format-number(TD[$xycont], '##0.0')"/></elep>
          </xsl:when>
          <xsl:otherwise>
            <elen><xsl:value-of select="format-number(TD[$xycont], '##0.0')"/></elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
<!-- Chi2 --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$chi2ref], '######0.0')"/></el>
        </td>
<!-- Shear --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$shear], '##0.00000')"/></el>
        </td>
<!-- Shear Position Angle--> 
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$sposangle], '##0.0')"/>&deg;</el>
        </td>
<!-- Zero-Point correction --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="contains(TD[$zpcorr],'e')">
           <el><xsl:value-of select="format-number('0', '#0.000')"/></el>
          </xsl:when>
          <xsl:otherwise>
           <el><xsl:value-of select="format-number(TD[$zpcorr], '#0.000')"/></el>
          </xsl:otherwise>
         </xsl:choose>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

<!-- ********************** XSL template for Fgroups ********************** -->
  <xsl:template name="FGroups">
   <xsl:variable name="name" select="count(FIELD[@name='Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="index" select="count(FIELD[@name='Index']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nfields" select="count(FIELD[@name='NFields']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="field_coord" select="count(FIELD[@name='Field_Coordinates']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pix_scale" select="count(FIELD[@name='Pixel_Scale']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="maxr" select="count(FIELD[@name='Max_Radius']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astrefcat" select="count(FIELD[@name='AstRef_Catalog']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astrefband" select="count(FIELD[@name='AstRef_Band']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astsigint" select="count(FIELD[@name='AstromSigma_Internal']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astcorrint" select="count(FIELD[@name='AstromCorr_Internal']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astchi2int" select="count(FIELD[@name='AstromChi2_Internal']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astndetint" select="count(FIELD[@name='AstromNDets_Internal']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astsigintH" select="count(FIELD[@name='AstromSigma_Internal_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astcorrintH" select="count(FIELD[@name='AstromCorr_Internal_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astchi2intH" select="count(FIELD[@name='AstromChi2_Internal_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astndetintH" select="count(FIELD[@name='AstromNDets_Internal_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astroffref" select="count(FIELD[@name='AstromOffset_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astsigref" select="count(FIELD[@name='AstromSigma_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astcorrref" select="count(FIELD[@name='AstromCorr_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astchi2ref" select="count(FIELD[@name='AstromChi2_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astndetref" select="count(FIELD[@name='AstromNDets_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astroffrefH" select="count(FIELD[@name='AstromOffset_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astsigrefH" select="count(FIELD[@name='AstromSigma_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astcorrrefH" select="count(FIELD[@name='AstromCorr_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astchi2refH" select="count(FIELD[@name='AstromChi2_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astndetrefH" select="count(FIELD[@name='AstromNDets_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nphot" select="PARAM[@name='NPhotInstru']"/>
   <xsl:variable name="photname" select="count(FIELD[@name='PhotInstru_Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="psigint" select="count(FIELD[@name='PhotSigma_Internal']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pchi2int" select="count(FIELD[@name='PhotChi2_Internal']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pndetint" select="count(FIELD[@name='PhotNDets_Internal']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="psigintH" select="count(FIELD[@name='PhotSigma_Internal_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pchi2intH" select="count(FIELD[@name='PhotChi2_Internal_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pndetintH" select="count(FIELD[@name='PhotNDets_Internal_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="psigref" select="count(FIELD[@name='PhotSigma_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pchi2ref" select="count(FIELD[@name='PhotChi2_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pndetref" select="count(FIELD[@name='PhotNDets_Reference']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="psigrefH" select="count(FIELD[@name='PhotSigma_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pchi2refH" select="count(FIELD[@name='PhotChi2_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pndetrefH" select="count(FIELD[@name='PhotNDets_Reference_HighSN']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fgplotflag" select="count(FIELD[@name='FgroupsPlot'])"/>
   <xsl:variable name="fgplot" select="count(FIELD[@name='FgroupsPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2plotflag" select="count(FIELD[@name='Chi2Plot'])"/>
   <xsl:variable name="chi2plot" select="count(FIELD[@name='Chi2Plot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ad1dplotflag" select="count(FIELD[@name='IntErr1DimPlot'])"/>
   <xsl:variable name="ad1dplot" select="count(FIELD[@name='IntErr1DimPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ad2dplotflag" select="count(FIELD[@name='IntErr2DimPlot'])"/>
   <xsl:variable name="ad2dplot" select="count(FIELD[@name='IntErr2DimPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="r1dplotflag" select="count(FIELD[@name='RefErr1DimPlot'])"/>
   <xsl:variable name="r1dplot" select="count(FIELD[@name='RefErr1DimPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="r2dplotflag" select="count(FIELD[@name='RefErr2DimPlot'])"/>
   <xsl:variable name="r2dplot" select="count(FIELD[@name='RefErr2DimPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pherplotflag" select="count(FIELD[@name='PhotErrPlot'])"/>
   <xsl:variable name="pherplot" select="count(FIELD[@name='PhotErrPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="phermagplotflag" select="count(FIELD[@name='PhotErrMagPlot'])"/>
   <xsl:variable name="phermagplot" select="count(FIELD[@name='PhotErrMagPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ZPplotflag" select="count(FIELD[@name='PhotZPPlot'])"/>
   <xsl:variable name="ZPplot" select="count(FIELD[@name='PhotZPPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ZP3plotflag" select="count(FIELD[@name='PhotZP3DPlot'])"/>
   <xsl:variable name="ZP3plot" select="count(FIELD[@name='PhotZP3DPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="csplotflag" select="count(FIELD[@name='ColShiftPlot'])"/>
   <xsl:variable name="csplot" select="count(FIELD[@name='ColShiftPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="mpplotflag" select="count(FIELD[@name='RefPropPlot'])"/>
   <xsl:variable name="mpplot" select="count(FIELD[@name='RefPropPlot']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('fgroups')">
     Groups Properties
    </BUTTON>
    <TABLE class="sortable" id="fgroups" BORDER="2" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Group Name</TH>
      <xsl:if test="$fgplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Group Plot</TH>
      </xsl:if>      
      <TH BGCOLOR="#FFEECC">Index</TH>
      <TH BGCOLOR="#FFEECC">Nfields</TH>
      <TH BGCOLOR="#FFEECC">Alpha</TH>
      <TH BGCOLOR="#FFEECC">Delta</TH>
      <TH BGCOLOR="#FFEECC">Pixel Scale</TH>
      <TH BGCOLOR="#FFEECC">Maximum Radius</TH>
      <TH BGCOLOR="#FFEECC">Reference Astrometric Catalog</TH>
      <TH BGCOLOR="#FFEECC">Reference Astrometric Band</TH>
      <xsl:if test="$chi2plotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">&chi;<sup>2</sup> Plot</TH>
      </xsl:if>      
      <TH BGCOLOR="#FFEECC">Internal Astrometric Sigma</TH>
      <TH BGCOLOR="#FFEECC">Internal Astrometric Correlation</TH>
      <TH BGCOLOR="#FFEECC">Internal &chi;<sup>2</sup></TH>
      <TH BGCOLOR="#FFEECC">Detections for internal statistics</TH>
      <TH BGCOLOR="#FFEECC">High S/N Internal Astrometric Sigma</TH>
      <TH BGCOLOR="#FFEECC">High S/N Internal Astrometric Correlation</TH>
      <TH BGCOLOR="#FFEECC">High S/N Internal &chi;<sup>2</sup></TH>
      <TH BGCOLOR="#FFEECC">High S/N Detections for Internal statistics</TH>
      <xsl:if test="$ad1dplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">1-D Internal Error Plot</TH>
      </xsl:if>      
      <xsl:if test="$ad2dplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">2-D Internal Error Plot</TH>
      </xsl:if>      
      <TH BGCOLOR="#FFEECC">Offset from Reference Catalog</TH>
      <TH BGCOLOR="#FFEECC">Sigma <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">Correlation <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">&chi;<sup>2</sup> <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">Detections for statistics <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N Offset from Reference Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N Sigma <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N Correlation <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N &chi;<sup>2</sup> <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N Detections for statistics <i>w.r.t.</i> Astrometric Ref Catalog</TH>
      <xsl:if test="$r1dplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">1-D External Error Plot</TH>
      </xsl:if>      
      <xsl:if test="$r2dplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">2-D External Error Plot</TH>
      </xsl:if>      
      <TH BGCOLOR="#FFEECC">Instrument Name</TH>
      <TH BGCOLOR="#FFEECC">Internal Photometric Sigma</TH>
      <TH BGCOLOR="#FFEECC">Internal &chi;<sup>2</sup></TH>
      <TH BGCOLOR="#FFEECC">Detections for internal statistics</TH>
      <TH BGCOLOR="#FFEECC">High S/N Internal Photometric Sigma</TH>
      <TH BGCOLOR="#FFEECC">High S/N Internal &chi;<sup>2</sup></TH>
      <TH BGCOLOR="#FFEECC">High S/N Detections for internal statistics</TH>
      <TH BGCOLOR="#FFEECC">Sigma <i>w.r.t.</i> Photometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">&chi;<sup>2</sup> <i>w.r.t.</i> Photometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">Detections for statistics <i>w.r.t.</i> Photometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N Sigma <i>w.r.t.</i> Photometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N &chi;<sup>2</sup> <i>w.r.t.</i> Photometric Ref Catalog</TH>
      <TH BGCOLOR="#FFEECC">High S/N Detections for statistics <i>w.r.t.</i> Photometric Ref Catalog</TH>
      <xsl:if test="$pherplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Internal Photometric Error Plot</TH>
      </xsl:if>      
      <xsl:if test="$phermagplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Internal Photometric Error vs. Magnitude Plot</TH>
      </xsl:if>      
      <xsl:if test="$ZPplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Zero Point Differences Plot</TH>
      </xsl:if>      
      <xsl:if test="$ZP3plotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">3D-Zero Point Differences Plot</TH>
      </xsl:if>      
      <xsl:if test="$csplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Color Shift Plot</TH>
      </xsl:if>      
      <xsl:if test="$mpplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Proper Motion Plot</TH>
      </xsl:if>      
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
        <xsl:if test="$fgplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$fgplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$fgplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$fgplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$index]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nfields]"/></el>
        </td>
<!-- Alpha -->
        <td align="center" BGCOLOR="#EEEEEE">
         <el>
          <xsl:variable name="alpha" select="substring-before(string(TD[$field_coord]), ' ')"/>
          <xsl:value-of
		select='concat(format-number(floor($alpha div 15.0), "00"),":",
		format-number(floor(($alpha * 4) mod 60.0), "00"),":",
		format-number(floor(($alpha * 240.0) mod 60.0), "00.00"))'/>
         </el>
        </td>
<!-- Delta -->
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:variable name="delta" select="substring-after(string(TD[$field_coord]), ' ')"/>
         <el>
          <xsl:choose>
           <xsl:when test="$delta &lt; 0.0">
            <xsl:value-of
		select='concat("-", format-number(floor(-$delta), "00"),":",
		format-number(floor((-$delta * 60) mod 60.0), "00"),":",
		format-number(floor((-$delta * 3600.0) mod 60.0), "00.0"))'/>
           </xsl:when>
           <xsl:otherwise>
            <xsl:value-of
		select='concat("+", format-number(floor($delta), "00"),":",
		format-number(floor(($delta * 60) mod 60.0), "00"),":",
		format-number(floor(($delta * 3600.0) mod 60.0), "00.0"))'/>
           </xsl:otherwise>
          </xsl:choose>
         </el>
        </td>
<!-- Pixel scale --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="pix1" select="number(substring-before(string(TD[$pix_scale]), ' '))"/>
         <xsl:variable name="pix2" select="number(substring-after(string(TD[$pix_scale]), ' '))"/>
         <el><xsl:value-of select="format-number(($pix1+$pix2) div 2.0, '##0.0000')"/>&asec;</el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el>
          <xsl:choose>
           <xsl:when test="TD[$maxr] &gt; 60.0">
            <xsl:value-of
             select="format-number(floor((TD[$maxr]) div 60.0), '##0.000')"/>&deg;
           </xsl:when>
           <xsl:otherwise>
            <xsl:value-of
             select="format-number(floor(TD[$maxr]), '#0.000')"/>&amin;
           </xsl:otherwise>
          </xsl:choose>
        </el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astrefcat]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astrefband]"/></el>
        </td>
        <xsl:if test="$chi2plotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$chi2plot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$chi2plot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$chi2plot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="sig1" select="number(substring-before(string(TD[$astsigint]), ' '))"/>
         <xsl:variable name="sig2" select="number(substring-after(string(TD[$astsigint]), ' '))"/>
         <el><xsl:value-of select='concat(format-number($sig1, "#0.0000"),"&asec; ",
                           format-number($sig2, "#0.0000"))'/>&asec;</el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astcorrint], '#0.00000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astchi2int], '####0.0')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astndetint]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="sigH1" select="number(substring-before(string(TD[$astsigintH]), ' '))"/>
         <xsl:variable name="sigH2" select="number(substring-after(string(TD[$astsigintH]), ' '))"/>
         <el><xsl:value-of select='concat(format-number($sigH1, "#0.0000"),"&asec; ",
                           format-number($sigH2, "#0.0000"))'/>&asec;</el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astcorrintH], '#0.00000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astchi2intH], '####0.0')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astndetintH]"/></el>
        </td>
        <xsl:if test="$ad1dplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$ad1dplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$ad1dplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$ad1dplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <xsl:if test="$ad2dplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$ad2dplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$ad2dplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$ad2dplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="off1" select="number(substring-before(string(TD[$astroffref]), ' '))"/>
         <xsl:variable name="off2" select="number(substring-after(string(TD[$astroffref]), ' '))"/>
         <el><xsl:value-of select='concat(format-number($off1, "#0.0000"),"&asec; ",
                           format-number($off2, "#0.0000"))'/>&asec;</el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="rsig1" select="number(substring-before(string(TD[$astsigref]), ' '))"/>
         <xsl:variable name="rsig2" select="number(substring-after(string(TD[$astsigref]), ' '))"/>
         <el><xsl:value-of select='concat(format-number($rsig1, "#0.000"),"&asec; ",
                           format-number($rsig2, "#0.000"))'/>&asec;</el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astcorrref], '#0.0000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astchi2ref], '####0.0')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astndetref]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="offH1" select="number(substring-before(string(TD[$astroffrefH]), ' '))"/>
         <xsl:variable name="offH2" select="number(substring-after(string(TD[$astroffrefH]), ' '))"/>
         <el><xsl:value-of select='concat(format-number($offH1, "#0.0000"),"&asec; ",
                           format-number($offH2, "#0.0000"))'/>&asec;</el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="rsigH1" select="number(substring-before(string(TD[$astsigrefH]), ' '))"/>
         <xsl:variable name="rsigH2" select="number(substring-after(string(TD[$astsigrefH]), ' '))"/>
         <el><xsl:value-of select='concat(format-number($rsigH1, "#0.000"),"&asec; ",
                           format-number($rsigH2, "#0.000"))'/>&asec;</el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astcorrrefH], '#0.0000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$astchi2refH], '####0.0')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astndetrefH]"/></el>
        </td>
        <xsl:if test="$r1dplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$r1dplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$r1dplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$r1dplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <xsl:if test="$r2dplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$r2dplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$r2dplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$r2dplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:value-of select="TD[$photname]"/>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$psigint]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pchi2int]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pndetint]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$psigintH]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pchi2intH]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pndetintH]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$psigref]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pchi2ref]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pndetref]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$psigrefH]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pchi2refH]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$pndetrefH]"/></el>
        </td>
        <xsl:if test="$pherplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$pherplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$pherplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$pherplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>      
        <xsl:if test="$phermagplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$phermagplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$phermagplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$phermagplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>      
        <xsl:if test="$ZPplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$ZPplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$ZPplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$ZPplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>      
        <xsl:if test="$ZP3plotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$ZP3plot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$ZP3plot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$ZP3plot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>      
        <xsl:if test="$csplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$csplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$csplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$csplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>      
        <xsl:if test="$mpplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$mpplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$mpplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$mpplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>      
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

<!-- ****************** XSL template for Astrometric Instruments ****************** -->
  <xsl:template name="Astr_Instr">
   <xsl:variable name="name" select="count(FIELD[@name='Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="index" select="count(FIELD[@name='Index']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nfields" select="count(FIELD[@name='NFields']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nkeys" select="count(FIELD[@name='NKeys']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="keys" select="count(FIELD[@name='Keys']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="distplotflag" select="count(FIELD[@name='DistPlot'])"/>
   <xsl:variable name="distplot" select="count(FIELD[@name='DistPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="refsysplotflag" select="count(FIELD[@name='RefSysPlot'])"/>
   <xsl:variable name="refsysplot" select="count(FIELD[@name='RefSysPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pixplotflag" select="count(FIELD[@name='PixErr1DimPlot'])"/>
   <xsl:variable name="pixplot" select="count(FIELD[@name='PixErr1DimPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="spixplotflag" select="count(FIELD[@name='SubPixErr1DimPlot'])"/>
   <xsl:variable name="spixplot" select="count(FIELD[@name='SubPixErr1DimPlot']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="shplotflag" select="count(FIELD[@name='ShearPlot'])"/>
   <xsl:variable name="shplot" select="count(FIELD[@name='ShearPlot']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('astro')">
     Astrometric Instruments
    </BUTTON>
    <TABLE class="sortable" id="astro" BORDER="2" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Name</TH>
      <TH BGCOLOR="#FFEECC">Index</TH>
      <TH BGCOLOR="#FFEECC">Nfields</TH>
      <TH BGCOLOR="#FFEECC">Number of Keywords</TH>
      <TH BGCOLOR="#FFEECC">Keywords</TH>
      <xsl:if test="$distplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Distortion Plot</TH>
      </xsl:if>      
      <xsl:if test="$refsysplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Systematics <i>w.r.t.</i> Reference Catalog Plot</TH>
      </xsl:if>      
      <xsl:if test="$pixplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Pixel Error Plot</TH>
      </xsl:if>      
      <xsl:if test="$spixplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Sub Pixel Error Plot</TH>
      </xsl:if>      
      <xsl:if test="$shplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Shear Plot</TH>
      </xsl:if>      
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$index]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nfields]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nkeys]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$keys]"/></el>
        </td>
        <xsl:if test="$distplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$distplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$distplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$distplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <xsl:if test="$refsysplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$refsysplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$refsysplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$refsysplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <xsl:if test="$pixplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$pixplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$pixplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$pixplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <xsl:if test="$spixplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$spixplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$spixplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$spixplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
        <xsl:if test="$shplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$shplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$shplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$shplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
       </tr>
      </xsl:for-each>
     </xsl:for-each>      
    </TABLE>
   </p>
 </xsl:template>

<!-- ****************** XSL template for Photometric Instruments ****************** -->
  <xsl:template name="Phot_Instr">
   <xsl:variable name="name" select="count(FIELD[@name='Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="index" select="count(FIELD[@name='Index']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nfields" select="count(FIELD[@name='NFields']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="magzerop" select="count(FIELD[@name='MagZeroPoint_Output']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nkeys" select="count(FIELD[@name='NKeys']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="keys" select="count(FIELD[@name='Keys']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('phot')">
     Photometric Instruments
    </BUTTON>
    <TABLE class="sortable" id="phot" BORDER="2" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Name</TH>
      <TH BGCOLOR="#FFEECC">Index</TH>
      <TH BGCOLOR="#FFEECC">Nfields</TH>
      <TH BGCOLOR="#FFEECC">Output ZP</TH>
      <TH BGCOLOR="#FFEECC">Number of Keywords</TH>
      <TH BGCOLOR="#FFEECC">Keywords</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$index]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nfields]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$magzerop]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nkeys]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$keys]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>      
    </TABLE>
   </p>
 </xsl:template>

<!-- ********************** XSL template for Config File ********************** -->
  <xsl:template name="Config">
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('config')">
     Configuration File: <xsl:value-of select="PARAM[@name='Prefs_Name']/@value"/>
    </BUTTON>
    <TABLE id="config" class="sortable" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Config Parameter</TH>
      <TH BGCOLOR="#FFEECC">Value</TH>
     </TR>
     <xsl:for-each select="PARAM[position()>2]">
      <tr BGCOLOR="#EEEEEE">
       <td><el><xsl:value-of select="@name"/></el></td>
       <td><el><xsl:value-of select="@value"/></el></td>
      </tr>
     </xsl:for-each>
    </TABLE>
   </p>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('commandline')">
     Command Line
    </BUTTON>
    <TABLE id="commandline" style="display: none">
     <TR>
      <TD BGCOLOR="#FFEECC" style="font-size: 80%;"><el>Command Line: <xsl:value-of select="PARAM[@name='Command_Line']/@value"/></el></TD>
     </TR>
    </TABLE>
   </p>
  </xsl:template>

<!-- ********************** XSL template for Warnings ********************** -->
  <xsl:template name="Warnings">
   <xsl:variable name="date" select="count(FIELD[@name='Date']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="time" select="count(FIELD[@name='Time']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="msg" select="count(FIELD[@name='Msg']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('warnings')">
     Warnings (limited to the last 1000)
    </BUTTON>
    <TABLE id="warnings" style="display: none">
     <TR style="font-size: 80%;">
      <TH BGCOLOR="#FFEECC">Date</TH>
      <TH BGCOLOR="#FFEECC">Time</TH>
      <TH BGCOLOR="#FFEECC">Message</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$msg]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

 <xsl:template name="Rest">
</xsl:template>
</xsl:stylesheet>
