from Bio import Entrez
import json
import time
import xml.etree.ElementTree as ET
import os

class Pubmed_API_langchain:
    def __init__(self, email):
        self.email = email
        Entrez.email = email
        self.papers_dict = {}  # Store papers in memory

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
        print(f"Title: {title}")
        
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

    def save_paper_to_memory(self, paper):
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
        
        # Store the paper data in memory
        self.papers_dict[paper['pmc_id']] = content
        
        print(f"Saved paper PMC{paper['pmc_id']} to memory")

    def query(self, query, max_results=10):
        print(f"PMC_API.query called with: {query}, {max_results}")

        id_list = self.search_pmc(query, max_results)
        papers = []
        for pmc_id in id_list:
            try:
                xml_content = self.fetch_full_text(pmc_id)
                paper_data = self.parse_full_text(xml_content)
                paper_data["pmc_id"] = pmc_id
                paper_data["citation_count"] = self.get_citation_count(pmc_id)
                papers.append(paper_data)
                
                # Save each paper to memory
                self.save_paper_to_memory(paper_data)
            except Exception as e:
                print(f"Error processing paper (PMC ID: {pmc_id}): {e}")
            time.sleep(1)  # Be nice to NCBI servers
        return papers

    def get_papers_for_llm(self):
        # Format papers for LLM agents
        formatted_papers = []
        for pmc_id, content in self.papers_dict.items():
            formatted_papers.append(f"Paper PMC{pmc_id}:\n\n{content}\n\n")
        return "\n".join(formatted_papers)

# Test the PMC_API
if __name__ == "__main__":
    pmc_api = Pubmed_API_langchain(email="sanazkazemi@hotmail.com")  # Replace with your email
    query = input("Enter a search query: ")
    max_results = 10

    print(f"Searching PMC for: '{query}'")
    papers = pmc_api.query(query, max_results)
    
    print(f"\nRetrieved and saved {len(papers)} papers in memory")

    # Example of how to get formatted papers for LLM agents
    llm_input = pmc_api.get_papers_for_llm()
    print("\nFormatted papers for LLM agents:")
    print(llm_input[:500] + "...")  # Print first 500 characters as an example

    # Here you would typically pass `llm_input` to your LLM agent
    # For example: llm_agent.process(llm_input)