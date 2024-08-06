from Bio import Entrez
import json
import time
import xml.etree.ElementTree as ET
import os
from paperqa import Docs, Answer, PromptCollection
from paperqa.prompts import select_paper_prompt, citation_prompt
from dotenv import load_dotenv
import tempfile
import asyncio

class PubMedProcessor:
    def __init__(self, email):
        load_dotenv()
        self.email = email
        Entrez.email = email
        self.openai_api_key = os.getenv("OPENAI_API_KEY")
        self.docs = None
        self.dictionary_for_llm = {}

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
        abstract = " ".join([get_text(p) for p in article_meta.findall(".//abstract//p")])
        authors = ", ".join([
            get_text(author.find(".//surname")) + " " + get_text(author.find(".//given-names"))
            for author in article_meta.findall(".//contrib[@contrib-type='author']")
        ])
        journal = get_text(root.find(".//journal-title"))
        year = get_text(article_meta.find(".//pub-date/year"))
        full_text = " ".join([get_text(p) for p in root.findall(".//body//p")])
        doi_elem = article_meta.find(".//article-id[@pub-id-type='doi']")
        doi = get_text(doi_elem) if doi_elem is not None else "No DOI available"
    
        return {
            "title": title or "No title available",
            "abstract": abstract,
            "authors": authors,
            "journal": journal or "Unknown Journal",
            "year": year or "N/A",
            "full_text": full_text,
            "doi": doi
        }

    def query_pubmed(self, query, max_results=10):
        id_list = self.search_pmc(query, max_results)
        papers = []
        for pmc_id in id_list:
            try:
                xml_content = self.fetch_full_text(pmc_id)
                paper_data = self.parse_full_text(xml_content)
                paper_data["pmc_id"] = pmc_id
                paper_data["citation_count"] = self.get_citation_count(pmc_id)
                papers.append(paper_data)
            except Exception as e:
                print(f"Error processing paper (PMC ID: {pmc_id}): {e}")
            time.sleep(1)  # Be nice to NCBI servers
        return papers

    def process_papers(self, papers):
        with tempfile.TemporaryDirectory() as temp_dir:
            prompts = PromptCollection()
            self.docs = Docs(prompts=prompts, llm="gpt-4o-mini")
            
            for i, paper in enumerate(papers):
                temp_file_path = os.path.join(temp_dir, f"temp_doc_{i}.txt")
                content = f"Title: {paper['title']}\n\n"
                content += f"Authors: {paper['authors']}\n\n"
                content += f"Journal: {paper['journal']}\n"
                content += f"Year: {paper['year']}\n"
                content += f"PMC ID: {paper['pmc_id']}\n"
                content += f"DOI: {paper['doi']}\n"
                content += f"Citation Count: {paper['citation_count']}\n\n"
                content += f"Abstract:\n{paper['abstract']}\n\n"
                content += f"Full Text:\n{paper['full_text']}\n"
                
                with open(temp_file_path, 'w', encoding='utf-8') as f:
                    f.write(content)
                
                try:
                    self.docs.add(temp_file_path)
                except Exception as e:
                    print(f"Error adding document to paperqa: {e}")

    def query_docs(self, query):
        if not self.docs:
            raise ValueError("Documents have not been processed. Call process_papers() first.")
        answer = self.docs.query(query)
        return answer

    def create_dictionary_for_llm(self, answer):
        evidence = self.docs.get_evidence(answer, k=10, max_sources=10, detailed_citations=True)
        self.dictionary_for_llm = {}

        for context in evidence.contexts:
            summary = f"summary: {context.context}"
            source_info = {
                "chunk_id": context.text.name,
                "full_citation": context.text.doc.citation,
            }
            chunk_info = {
                "original_text": context.text.text,
                "source": source_info,
                "relevance_score": context.score
            }
            self.dictionary_for_llm[summary] = chunk_info

        return self.dictionary_for_llm

    async def full_process(self, pubmed_query, doc_query, max_results=10):
        papers = self.query_pubmed(pubmed_query, max_results)
        await asyncio.to_thread(self.process_papers, papers)
        answer = self.query_docs(doc_query)
        self.dictionary_for_llm = self.create_dictionary_for_llm(answer)

        return self.create_dictionary_for_llm(answer)

    # def get_dictionary(self): 
    #     return self.dictionary_for_llm

# Usage example:
# if __name__ == "__main__":
#     async def main():
#         processor = PubMedProcessor(email="your_email@example.com")
#         pubmed_query = "dopamine 1 receptor Parkinson's disease"
#         doc_query = "relevance of dopamine 1 receptor to Parkinson's disease"
#         result = await processor.full_process(pubmed_query, doc_query, max_results=5)
#         print(json.dumps(result, indent=2))

#     asyncio.run(main())

        
