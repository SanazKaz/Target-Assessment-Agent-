from Bio import Entrez
import json
import time
import xml.etree.ElementTree as ET
import os

class Pubmed_API:
    def __init__(self, email):
        self.email = email
        Entrez.email = email

    def search_pmc(self, query, max_results=10):
        handle = Entrez.esearch(db="pmc", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]

    def get_citation_count(self, pmc_id):
        handle = Entrez.elink(dbfrom="pmc", db="pmc", linkname="pmc_pmc_citedby", id=pmc_id)
        result = Entrez.read(handle)
        handle.close()
        if result[0]["LinkSetDb"]:
            return len(result[0]["LinkSetDb"][0]["Link"])
        return 0

    def fetch_full_text(self, pmc_id):
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="xml")
        xml_content = handle.read()
        handle.close()
        return xml_content

    def parse_full_text(self, xml_content):
        root = ET.fromstring(xml_content)
        article_meta = root.find(".//article-meta")
    
        def get_text(element):
            if element is None:
                return ""
            return ''.join(element.itertext()).strip()

        title = get_text(article_meta.find(".//article-title"))
        if not title:
            title = "No title available"
        
        abstract = ""
        abstract_elem = article_meta.find(".//abstract")
        if abstract_elem is not None:
            abstract = " ".join([get_text(p) for p in abstract_elem.findall(".//p")])
        
        authors = ", ".join([
            get_text(author.find(".//surname")) + " " + get_text(author.find(".//given-names"))
            for author in article_meta.findall(".//contrib[@contrib-type='author']")
        ])
        
        journal = get_text(root.find(".//journal-title"))
        if not journal:
            journal = "Unknown Journal"
        
        year = get_text(article_meta.find(".//pub-date/year"))
        if not year:
            year = "N/A"
        
        full_text = " ".join([get_text(p) for p in root.findall(".//body//p")])
        
        doi_elem = article_meta.find(".//article-id[@pub-id-type='doi']")
        doi = get_text(doi_elem) if doi_elem is not None else "No DOI available"
    
        return {
            "title": title,
            "abstract": abstract,
            "authors": authors,
            "journal": journal,
            "year": year,
            "full_text": full_text,
            "doi": doi
        }

    def save_paper_to_txt(self, paper, output_dir="pmc_papers"):
        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Create a filename using the PMC ID
        filename = os.path.join(output_dir, f"PMC{paper['pmc_id']}.txt")
        
        # Format the paper data as text
        content = f"Title: {paper['title']}\n\n"
        content += f"Authors: {paper['authors']}\n\n"
        content += f"Journal: {paper['journal']}\n"
        content += f"Year: {paper['year']}\n"
        content += f"PMC ID: {paper['pmc_id']}\n"
        content += f"DOI: {paper['doi']}\n"
        content += f"Citation Count: {paper['citation_count']}\n\n"
        content += f"Abstract:\n{paper['abstract']}\n\n"
        content += f"Full Text:\n{paper['full_text']}\n"
        
        # Save the paper data to a text file
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f"Saved paper to {filename}")

    def query(self, query, max_results=10):
        print(f"PMC_API.query called with: {query}, {max_results}")

        id_list = self.search_pmc(query, max_results)  # Get results sorted by relevance
        papers = []
        for pmc_id in id_list:
            try:
                xml_content = self.fetch_full_text(pmc_id)
                paper_data = self.parse_full_text(xml_content)
                paper_data["pmc_id"] = pmc_id
                paper_data["citation_count"] = self.get_citation_count(pmc_id)
                papers.append(paper_data)
                
                # Save each paper to its own text file
                self.save_paper_to_txt(paper_data)
            except Exception as e:
                print(f"Error processing paper (PMC ID: {pmc_id}): {e}")
            time.sleep(1)  # Be nice to NCBI servers
        return papers

# Test the PMC_API
if __name__ == "__main__":
    pmc_api = Pubmed_API(email="sanazkazemi@hotmail.com")  # Replace with your email
    query = input("Enter a search query: ")
    max_results = 10

    print(f"Searching PMC for: '{query}'")
    papers = pmc_api.query(query, max_results)
    
    print(f"\nRetrieved and saved {len(papers)} papers:")
    for paper in papers:
        print(f"\nTitle: {paper['title']}")
        print(f"Authors: {paper['authors']}")
        print(f"Journal: {paper['journal']}, {paper['year']}")
        print(f"PMC ID: {paper['pmc_id']}")
        print(f"DOI: {paper['doi']}")
        print(f"Citation Count: {paper['citation_count']}")
        print(f"Abstract: {paper['abstract'][:200]}...")
        print(f"Full Text: {paper['full_text'][:200]}...")

    print("\nAll papers have been saved to individual text files in the 'pmc_papers' directory.")