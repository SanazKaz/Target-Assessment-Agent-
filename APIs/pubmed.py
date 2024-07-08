from Bio import Entrez
import time

class PubMedAPI:
    def __init__(self, email):
        self.email = email
        Entrez.email = email

    def search_pubmed(self, query, max_results=10):
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]

    def fetch_details(self, id_list):
        ids = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        results = Entrez.read(handle)
        handle.close()
        return results

    def query(self, query, max_results=10):
        id_list = self.search_pubmed(query, max_results)
        papers = []
        for paper in self.fetch_details(id_list)["PubmedArticle"]:
            title = paper["MedlineCitation"]["Article"]["ArticleTitle"]
            try:
                abstract = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
            except KeyError:
                abstract = "No abstract available"
            authors = ", ".join([author["LastName"] + " " + author["Initials"] 
                                 for author in paper["MedlineCitation"]["Article"]["AuthorList"]])
            journal = paper["MedlineCitation"]["Article"]["Journal"]["Title"]
            year = paper["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
            
            # Extract DOI
            doi = ""
            for id in paper["PubmedData"]["ArticleIdList"]:
                if id.attributes["IdType"] == "doi":
                    doi = str(id)
                    break
            
            papers.append({
                "title": title,
                "abstract": abstract,
                "authors": authors,
                "journal": journal,
                "year": year,
                "doi": doi
            })
            time.sleep(1)  # Be nice to NCBI servers
        return papers

    def format_results(self, papers):
        return "\n\n".join([
            f"Title: {paper['title']}\n"
            f"Authors: {paper['authors']}\n"
            f"Journal: {paper['journal']}, {paper['year']}\n"
            f"DOI: {paper['doi']}\n"
            f"Abstract: {paper['abstract'][:200]}..."
            for paper in papers
        ])

# Test the PubMedAPI
if __name__ == "__main__":
    pubmed_api = PubMedAPI(email="your_email@example.com")  # Replace with your email
    query = "ACE inhibitors hypertension"
    max_results = 5

    print(f"Searching PubMed for: '{query}'")
    papers = pubmed_api.query(query, max_results)
    
    print(f"\nRetrieved {len(papers)} papers:")
    for paper in papers:
        print(f"\nTitle: {paper['title']}")
        print(f"Authors: {paper['authors']}")
        print(f"Journal: {paper['journal']}, {paper['year']}")
        print(f"DOI: {paper['doi']}")
        print(f"Abstract: {paper['abstract'][:200]}...")

    print("\n\nFormatted results:")
    print(pubmed_api.format_results(papers))