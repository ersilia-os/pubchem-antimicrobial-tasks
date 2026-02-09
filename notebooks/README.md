# About PubChem 💊
### What is PubChem?
[PubChem](https://pubchem.ncbi.nlm.nih.gov/) is the world's largest collection of freely accessible chemical information. Search chemicals by name, molecular formula, structure, and other identifiers. Find chemical and physical properties, biological activities, safety and toxicity information, patents, literature citations and more.
### PubChem Bioassays 🧪
**PubChem BioAssay** contains small-molecule and RNAi screening data along with associated annotation information from contributing organizations. 

BioAssay contains a collection of **bioactivity** and **toxicity** data that has greatly supported research in fields such as medicinal chemistry, drug discovery, pharmaceutical genomics and informatics research.

Each **test result** provided points to an entry in the PubChem Substance database using a unique identifier for that molecule's record. If the tested substance has a known structure it is displayed. Each test result must provide an expert opinion on the overall activity of the substance in this experiment, which is typically **active** or **inactive**. 

To access a PubChem BioAssay record by its **numeric identifier (AID)**, e.g. https://pubchem.ncbi.nlm.nih.gov/bioassay/1. Or search in PubChem with some query (e.g. assay description, compound name, and protein/gene target) and select the BioAssays tab. 

For each BioAssay record, **bioactivity data** together with chemical structures (in isomeric **SMILES** format) are available for download, which might facilitate the studies of structure-activity analysis without the extra need to obtain chemical structures from PubChem CIDs.

## Downloading PubChem Bioassay Data ⬇️
[PubChem BioAssay](https://pubchem.ncbi.nlm.nih.gov/docs/bioassays#section=Programmatic-Access-Examples) data are available for **bulk download** on the PubChem FTP site (https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/).

[Programmatic services](https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access) can also be used to access PubChem BioAssay data ([PUG REST Tutorial](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial)).

### How PUG REST Works
The fundamental unit upon which PUG REST is built is the **PubChem identifier**, which comes in three flavors – **SID for substances**, **CID for compounds**, and **AID for assays**. 

The conceptual framework of this service, that uses these identifiers, is the three-part request: 

1) **input** – that is, what identifiers are we talking about; 
2) **operation** – what to do with those identifiers;
3) **output** – what information should be returned. 

The beauty of this design is that each of these three parts of the request is (mostly) independent, allowing a combinatorial expansion of the things you can do in a single request. Meaning that, for example, any form of input that specifies some group of CIDs can be combined with any operation that deals with CIDs, and any output format that’s relevant to the chosen operation. So instead of a list of separate narrowly defined service requests that are supported, you can combine these building blocks in many ways to create customized requests.

Example: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/vioxx/property/InChI/TXT

#### 1. Input: Design of the URL
PUG REST is entirely based on HTTP (or HTTPS) requests, and most of the details of the request are encoded directly in the URL path.

Analysing last example:

| prolog | input | operation | output |
|--------|--------|-----------|---------|
| https://pubchem.ncbi.nlm.nih.gov/rest/pug | compound/name/vioxx | property/InChI | TXT |

Taking each section individually, first we have the prolog – the HTTP address of the service itself, which is common to all PUG REST requests. The next part is the input, which in this case says “I want to look in the PubChem Compound database for records that match the name ‘vioxx’.” Note that there some subtleties here, in that the name must already be present in the PubChem database, and that a name may refer to multiple CIDs. But the underlying principle is that we are specifying a set of CIDs based on a name; at the time of writing, there is only one CID with this name. The next section is the operation, in this case “I want to retrieve the InChI property for this CID.” And finally the output format specification, “I want to get back plain text.”

#### 2. Output: What You Get Back
The results of most operations can be expressed in a variety of data formats:
| Output Format | Description |
|---------------|-------------|
| XML   | standard XML, for which a schema is available |
| JSON  | JSON, JavaScript Object Notation |
| JSONP | JSONP, like JSON but wrapped in a callback function |
| ASNB  | standard binary ASN.1, NCBI’s native format in many cases |
| ASNT  | NCBI’s human-readable text flavor of ASN.1 |
| SDF   | chemical structure data |
| CSV   | comma-separated values, spreadsheet compatible |
| PNG   | standard PNG image data |
| TXT   | plain text |

### Access to PubChem BioAssays (AID)
An assay is composed of two general parts: the **description** of the assay as a whole, including authorship, general description, protocol, and definitions of the data readout columns; and then the **data** section with all the actual test result values.

1) **Assay Description**

    To get just the description section via PUG REST, use a request like:

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/504526/description/XML

    There is also a simplified summary format that does not have the full complexity of the original description as above, and includes some information on targets, active and inactive SID and CID counts, etc. For example:

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/1000/summary/JSON

2) **Assay Data**

    BioAssay data may be conceptualized as a large table where the **columns** are the **readouts** (enumerated in the description section), and the **rows** are the **individual substances tested** and their results for each column. So, retrieving an entire assay record involves the primary AID – the identifier for the assay itself – and a list of SIDs. 
    
    If you want all the data rows of an assay, you can use a simple request like this one, which will return a CSV table of results. Note that full-data retrieval is the default operation for assays.

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/504526/CSV

    However, as some assays have many thousands of SID rows, there is a limit, currently 10,000, on the number of rows that can be retrieved in a single request. If you are interested in only a subset of the total data rows, you can use an optional argument to the PUG REST request to limit the output to just those SIDs (and note that with XML/ASN output you get the description as well when doing data retrieval). 
    
    There are other ways to input the SID list, such as in the HTTP POST body or via a list key; see below for more detail on lists stored on the server.

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/504526/XML?sid=104169547,109967232

3) Targets

    When the target of a BioAssay is known, it can be retrieved either as a sequence or gene, including identifiers in NCBI’s respective databases for these: https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/490,1000/targets/ProteinGI,ProteinName,GeneID,GeneSymbol/XML

4) Activity Name

    BioAssays may be selected by the name of the primary activity column, for example to get all the AIDs that are measuring an EC50:

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/activity/EC50/aids/JSON

### Access to PubChem Taxonomies (Taxonomy ID)
1) **Taxonomy Summary**

    This operation returns a summary of taxonomy: TaxonomyID, ScientificName, CommonName, Rank, RankedLineage, and a list of Synonyms. Valid output formats are XML, JSON(P), and ASNT/B. For example:

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/taxonomy/taxid/9606,10090,10116/summary/JSON

2) **Assays and Bioactivities**

    The following operation returns a **list of compounds** involved in a given taxonomy. Valid output formats are XML, JSON(P), ASNT/B, and TXT.

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/taxonomy/taxid/2697049/aids/TXT

    There is **no operation available to directly retrieve the bioactivity data associated with a given taxonomy**, as often the data volume is huge. However, one can first **get the list of AIDs** using the above link, and then aggregate the concise bioactivity data from each AID, e.g.:

    https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/1409578/concise/JSON

## PubChem Taxonomies 🪲
You can access a PubChem taxonomy via a URL of this form:

https://pubchem.ncbi.nlm.nih.gov/taxonomy/TAXON

where TAXON is the name of a taxon.  It can be a common name (e.g., rat or Norway rat) or a scientific name (Rattus norvegicus).