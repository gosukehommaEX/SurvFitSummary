"""Write report_template.docx, the Word document whose styles
write_survival_report() uses by default.

The document has no content (one empty paragraph); only its page set-up and
styles are used: A4 page with 2.5 cm margins, Arial 10 pt text,
numbered headings (1., 1.1., 1.1.1.) in dark red with a rule under
Heading 1, an unnumbered Heading 4, bold 9 pt captions, a List Bullet
style for bulleted lists, and tables with
grey grid lines and a dark red header row with white bold text. Table text
is 9 pt so that the wide tables fit the page.

Run from this folder with Python 3 (standard library only):
    python3 make_report_template.py
"""
import sys
import zipfile

W = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'
R = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
HEADER = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'

CONTENT_TYPES = HEADER + """<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Default Extension="png" ContentType="image/png"/>
<Default Extension="jpeg" ContentType="image/jpeg"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
<Override PartName="/word/numbering.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.numbering+xml"/>
<Override PartName="/word/settings.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.settings+xml"/>
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>"""

PACKAGE_RELS = HEADER + """<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>"""

DOCUMENT_RELS = HEADER + """<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/numbering" Target="numbering.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/settings" Target="settings.xml"/>
</Relationships>"""

CORE = HEADER + """<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<dc:title>Survival Analysis Report</dc:title>
<dc:creator>SurvFitSummary</dc:creator>
<cp:lastModifiedBy>SurvFitSummary</cp:lastModifiedBy>
<cp:revision>1</cp:revision>
<dcterms:created xsi:type="dcterms:W3CDTF">2026-11-13T00:00:00Z</dcterms:created>
<dcterms:modified xsi:type="dcterms:W3CDTF">2026-11-13T00:00:00Z</dcterms:modified>
</cp:coreProperties>"""

APP = HEADER + """<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
<Application>SurvFitSummary</Application>
</Properties>"""

SETTINGS = HEADER + f"""<w:settings {W}>
<w:defaultTabStop w:val="720"/>
<w:characterSpacingControl w:val="doNotCompress"/>
<w:compat><w:compatSetting w:name="compatibilityMode" w:uri="http://schemas.microsoft.com/office/word" w:val="15"/></w:compat>
</w:settings>"""

# A4 (11906 x 16838 twips) with 2.5 cm (1417 twips) margins
DOCUMENT = HEADER + f"""<w:document {W} {R}>
<w:body>
<w:p/>
<w:sectPr>
<w:pgSz w:w="11906" w:h="16838"/>
<w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" w:header="708" w:footer="708" w:gutter="0"/>
<w:cols w:space="708"/>
</w:sectPr>
</w:body>
</w:document>"""

# Headings 1 to 3 numbered 1., 1.1., 1.1.1. (numId 1), and the bullets of
# the style List Bullet (numId 2)
NUMBERING = HEADER + f"""<w:numbering {W}>
<w:abstractNum w:abstractNumId="0">
<w:multiLevelType w:val="multilevel"/>
<w:lvl w:ilvl="0"><w:start w:val="1"/><w:numFmt w:val="decimal"/><w:pStyle w:val="Heading1"/><w:lvlText w:val="%1."/><w:lvlJc w:val="left"/><w:pPr><w:ind w:left="432" w:hanging="432"/></w:pPr></w:lvl>
<w:lvl w:ilvl="1"><w:start w:val="1"/><w:numFmt w:val="decimal"/><w:pStyle w:val="Heading2"/><w:lvlText w:val="%1.%2."/><w:lvlJc w:val="left"/><w:pPr><w:ind w:left="576" w:hanging="576"/></w:pPr></w:lvl>
<w:lvl w:ilvl="2"><w:start w:val="1"/><w:numFmt w:val="decimal"/><w:pStyle w:val="Heading3"/><w:lvlText w:val="%1.%2.%3."/><w:lvlJc w:val="left"/><w:pPr><w:ind w:left="720" w:hanging="720"/></w:pPr></w:lvl>
</w:abstractNum>
<w:abstractNum w:abstractNumId="1">
<w:multiLevelType w:val="singleLevel"/>
<w:lvl w:ilvl="0"><w:start w:val="1"/><w:numFmt w:val="bullet"/><w:pStyle w:val="ListBullet"/><w:lvlText w:val="\u2022"/><w:lvlJc w:val="left"/><w:pPr><w:ind w:left="360" w:hanging="360"/></w:pPr><w:rPr><w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:cs="Arial"/></w:rPr></w:lvl>
</w:abstractNum>
<w:num w:numId="1"><w:abstractNumId w:val="0"/></w:num>
<w:num w:numId="2"><w:abstractNumId w:val="1"/></w:num>
</w:numbering>"""


def paragraph_style(style_id, name, ppr="", rpr="", based_on="Normal",
                    next_style="Normal", custom=False):
    custom_attr = ' w:customStyle="1"' if custom else ""
    return (f'<w:style w:type="paragraph"{custom_attr} w:styleId="{style_id}">'
            f'<w:name w:val="{name}"/><w:basedOn w:val="{based_on}"/>'
            f'<w:next w:val="{next_style}"/><w:qFormat/>'
            f'<w:pPr>{ppr}</w:pPr><w:rPr>{rpr}</w:rPr></w:style>')


def heading(level, size, color, before, after, numbered, border=""):
    num = (f'<w:numPr><w:ilvl w:val="{level - 1}"/><w:numId w:val="1"/></w:numPr>'
           if numbered else "")
    return paragraph_style(
        f"Heading{level}", f"heading {level}",
        ppr=(f'<w:keepNext/><w:keepLines/>{num}{border}'
             f'<w:spacing w:before="{before}" w:after="{after}"/>'
             f'<w:outlineLvl w:val="{level - 1}"/>'),
        rpr=f'<w:b/><w:bCs/><w:color w:val="{color}"/><w:sz w:val="{size}"/><w:szCs w:val="{size}"/>'
    )


GRID = "".join(f'<w:{side} w:val="single" w:sz="4" w:space="0" w:color="7F7F7F"/>'
               for side in ["top", "left", "bottom", "right", "insideH", "insideV"])

STYLES = HEADER + f"""<w:styles {W}>
<w:docDefaults>
<w:rPrDefault><w:rPr><w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:eastAsia="Arial" w:cs="Arial"/><w:sz w:val="20"/><w:szCs w:val="20"/><w:lang w:val="en-GB" w:eastAsia="en-US" w:bidi="ar-SA"/></w:rPr></w:rPrDefault>
<w:pPrDefault><w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr></w:pPrDefault>
</w:docDefaults>
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/><w:pPr><w:spacing w:before="120"/></w:pPr><w:rPr><w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:cs="Arial"/><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr></w:style>
<w:style w:type="character" w:default="1" w:styleId="DefaultParagraphFont"><w:name w:val="Default Paragraph Font"/><w:uiPriority w:val="1"/><w:semiHidden/><w:unhideWhenUsed/></w:style>
<w:style w:type="table" w:default="1" w:styleId="TableNormal"><w:name w:val="Normal Table"/><w:uiPriority w:val="99"/><w:semiHidden/><w:unhideWhenUsed/><w:tblPr><w:tblInd w:w="0" w:type="dxa"/><w:tblCellMar><w:top w:w="0" w:type="dxa"/><w:left w:w="108" w:type="dxa"/><w:bottom w:w="0" w:type="dxa"/><w:right w:w="108" w:type="dxa"/></w:tblCellMar></w:tblPr></w:style>
<w:style w:type="numbering" w:default="1" w:styleId="NoList"><w:name w:val="No List"/><w:uiPriority w:val="99"/><w:semiHidden/><w:unhideWhenUsed/></w:style>
{paragraph_style("Title", "Title", ppr='<w:spacing w:before="0" w:after="120"/>', rpr='<w:b/><w:bCs/><w:sz w:val="44"/><w:szCs w:val="44"/>')}
{paragraph_style("Subtitle", "Subtitle", ppr='<w:spacing w:before="0" w:after="240"/>', rpr='<w:i/><w:iCs/><w:color w:val="7F7F7F"/><w:sz w:val="22"/><w:szCs w:val="22"/>')}
{heading(1, 32, "C00000", 360, 120, True, border='<w:pBdr><w:bottom w:val="single" w:sz="4" w:space="1" w:color="auto"/></w:pBdr>')}
{heading(2, 26, "C00000", 240, 80, True)}
{heading(3, 22, "C00000", 200, 60, True)}
{heading(4, 20, "800000", 160, 40, False)}
{paragraph_style("Caption", "caption", ppr='<w:spacing w:before="200" w:after="60"/>', rpr='<w:b/><w:bCs/><w:sz w:val="18"/><w:szCs w:val="18"/>')}
{paragraph_style("TableCaption", "Table Caption", ppr='<w:keepNext/><w:keepLines/>', based_on="Caption", custom=True)}
{paragraph_style("Figure", "Figure", ppr='<w:keepNext/><w:spacing w:before="200"/>', custom=True)}
{paragraph_style("TableText", "Table Text", ppr='<w:spacing w:before="20" w:after="20"/>', rpr='<w:sz w:val="18"/><w:szCs w:val="18"/>', custom=True)}
{paragraph_style("ListBullet", "List Bullet", ppr='<w:numPr><w:ilvl w:val="0"/><w:numId w:val="2"/></w:numPr><w:spacing w:before="40"/><w:ind w:left="360" w:hanging="360"/>')}
<w:style w:type="table" w:customStyle="1" w:styleId="ReportTable"><w:name w:val="Report Table"/><w:basedOn w:val="TableNormal"/><w:uiPriority w:val="59"/><w:tblPr><w:tblBorders>{GRID}</w:tblBorders><w:tblCellMar><w:left w:w="57" w:type="dxa"/><w:right w:w="57" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblStylePr w:type="firstRow"><w:rPr><w:b/><w:bCs/><w:color w:val="FFFFFF"/></w:rPr><w:tcPr><w:shd w:val="clear" w:color="auto" w:fill="800000"/></w:tcPr></w:tblStylePr></w:style>
</w:styles>"""


def main(path="report_template.docx"):
    parts = [
        ("[Content_Types].xml", CONTENT_TYPES),
        ("_rels/.rels", PACKAGE_RELS),
        ("docProps/core.xml", CORE),
        ("docProps/app.xml", APP),
        ("word/_rels/document.xml.rels", DOCUMENT_RELS),
        ("word/document.xml", DOCUMENT),
        ("word/styles.xml", STYLES),
        ("word/numbering.xml", NUMBERING),
        ("word/settings.xml", SETTINGS),
    ]
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in parts:
            info = zipfile.ZipInfo(name, date_time=(2026, 11, 13, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, text.encode("utf-8"))


if __name__ == "__main__":
    main(*sys.argv[1:])
