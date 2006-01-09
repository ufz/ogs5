<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:ogs="http://www.opengeosys.net" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="/">
  <html>
  <body>
    <h4 style="font-family:verdana">OpenGeoSys - Geometrical data file</h4>
    <h2><xsl:value-of select="OpenGeoSysGLI/name" /></h2>
    
    <table cellpadding="30">
    
    <tr><td valign="top" align="left" width="33%">
    <p style="font-family:verdana"><strong>Points:</strong></p>
    <ul>
      <xsl:for-each select="OpenGeoSysGLI/points/point">
        <li>Point <xsl:value-of select="@id"/> <div style="font-family:arial;font-size:small"> ( <xsl:value-of select="@x"/>, <xsl:value-of select="@y"/> <xsl:if test="@z">, <xsl:value-of select="@z"/> </xsl:if> ) <p> </p></div></li>
      </xsl:for-each>
    </ul>
    
    </td><td valign="top" align="left" width="33%">
    <p style="font-family:verdana"><strong>Polylines:</strong></p>
    
    <ul>
      <xsl:for-each select="OpenGeoSysGLI/polylines/polyline">
        <li>Polyline <xsl:value-of select="@id"/>:<br />
        <div style="font-family:arial;font-size:small">
        <xsl:for-each select="./pnt">
          <xsl:apply-templates/> -
        </xsl:for-each>
        </div>
        <p> </p>
        </li>
      </xsl:for-each>
    </ul>
    
    </td><td valign="top" align="left" width="33%">
    <p style="font-family:verdana"><strong>Surfaces:</strong></p>
    
    <ul>
      <xsl:for-each select="OpenGeoSysGLI/surfaces/surface">
        <li>Surface <xsl:value-of select="@id"/>:<br />
<!--          
            (<xsl:choose>
              <xsl:when test="@matgroup"> 
                Material Group <xsl:value-of select="@matgroup"/> 
              </xsl:when>
              <xsl:otherwise>
                <em>no material group</em>
              </xsl:otherwise>
            </xsl:choose>
            , 
            <xsl:choose>
              <xsl:when test="@epsilon">
                Epsilion=<xsl:value-of select="@epsilon"/> 
              </xsl:when>
              <xsl:otherwise>
                <em>no epsilon value</em>
              </xsl:otherwise>
            </xsl:choose>
          ) 
-->
	  <ul>
            <xsl:for-each select="./element">
	      <li style="font-family:arial;font-size:small">Triangle( <xsl:value-of select="@p1"/>, <xsl:value-of select="@p2"/>, <xsl:value-of select="@p3"/> )</li>
	    </xsl:for-each>
	  </ul>
          <br />

        </li>
      </xsl:for-each>
    </ul>
    
    </td></tr>
    </table>
  </body>
  </html>
  
</xsl:template>

</xsl:stylesheet>