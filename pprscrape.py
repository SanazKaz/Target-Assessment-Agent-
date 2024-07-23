import requests
import xml.etree.ElementTree as ET
import time
import json
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

def requests_retry_session(
    retries=3,
    backoff_factor=0.3,
    status_forcelist=(500, 502, 504),
    session=None,
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def search_pmc(query, retmax=10):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pmc",
        "term": query,
        "retmax": retmax,
        "retmode": "json"
    }
    try:
        response = requests_retry_session().get(base_url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        return data['esearchresult']['idlist']
    except Exception as e:
        print(f"An error occurred while searching: {e}")
        return []

def get_pmc_fulltext(pmcid):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pmc",
        "id": pmcid,
        "retmode": "xml"
    }
    try:
        response = requests_retry_session().get(base_url, params=params, timeout=10)
        response.raise_for_status()
        root = ET.fromstring(response.content)
        
        full_text = ' '.join(root.itertext())
        title = root.find(".//article-title").text if root.find(".//article-title") is not None else "No title found"
        abstract = ' '.join(root.find(".//abstract").itertext()) if root.find(".//abstract") is not None else "No abstract found"
        
        return {
            "pmcid": pmcid,
            "title": title,
            "abstract": abstract,
            "full_text": full_text
        }
    except Exception as e:
        print(f"An error occurred while retrieving full text for PMC{pmcid}: {e}")
        return None

query = "machine learning AND covid-19"  # Replace with your search query
pmcids = search_pmc(query, retmax=5)  # Get 5 results

articles = []

for pmcid in pmcids:
    print(f"Retrieving article {pmcid}")
    article_data = get_pmc_fulltext(pmcid)
    if article_data:
        articles.append(article_data)
        print(f"Retrieved data for PMC{pmcid}")
    else:
        print(f"Failed to retrieve data for PMC{pmcid}")
    
    # Longer delay between requests
    time.sleep(3)

# Save to JSON file
output_file = "pmc_articles.json"
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(articles, f, ensure_ascii=False, indent=4)

print(f"Saved {len(articles)} articles to {output_file}")