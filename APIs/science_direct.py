import requests
import json
from dotenv import load_dotenv
import os

load_dotenv()

INST_token = os.getenv("SCIENCE_DIRECT_INST_KEY")
API_KEY = os.getenv('SCIENCE_DIRECT_API_KEY')

SEARCH_URL = 'https://api.elsevier.com/content/search/sciencedirect'
ARTICLE_URL = 'https://api.elsevier.com/content/article/doi/'

def search_articles(query):
    params = {
        'query': query,
    }
    headers = {
        'X-ELS-APIKey': API_KEY,
        'X-ELS-Insttoken': INST_token,
        'Accept': 'application/json'
    }
    
    response = requests.get(SEARCH_URL, params=params, headers=headers)
    
    if response.status_code == 200:
        print(f"Search status code: {response.status_code}")
        return response.json()
    else:
        print(f"Search error: {response.status_code}")
        print(response.text)
        return None

def retrieve_article(doi):
    headers = {
        'X-ELS-APIKey': API_KEY,
        'X-ELS-Insttoken': INST_token,
        'Accept': 'application/json'
    }
    
    url = f"{ARTICLE_URL}{doi}"
    
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        print(f"Article retrieval status code: {response.status_code}")
        return response.json()
    else:
        print(f"Article retrieval error: {response.status_code}")
        print(response.text)
        return None

def display_search_results(results):
    if results and 'search-results' in results:
        entries = results['search-results'].get('entry', [])
        print(f"Found {len(entries)} articles:")
        for i, article in enumerate(entries, 1):
            title = article.get('dc:title', 'No title available')
            doi = article.get('prism:doi', 'No DOI available')
            print(f"{i}. Title: {title}")
            print(f"   DOI: {doi}")
            print()
        return entries
    else:
        print("No results found or unexpected response format.")
        return []

def display_article(article):
    if article and 'full-text-retrieval-response' in article:
        data = article['full-text-retrieval-response']
        title = data.get('coredata', {}).get('dc:title', 'No title available')
        authors = data.get('coredata', {}).get('dc:creator', 'No authors available')
        abstract = data.get('coredata', {}).get('dc:description', 'No abstract available')
        doi = data.get('coredata', {}).get('prism:doi', 'No DOI available')
        
        print(f"Title: {title}")
        print(f"Authors: {authors}")
        print(f"Abstract: {abstract}")
        print(f"doi:{doi}")
        
        # You can add more fields as needed
    else:
        print("Unexpected response format or no article data found.")

def main():
    search_query = input("Enter your search query: ")
    results = search_articles(search_query)
    
    if results:
        entries = display_search_results(results)
        if entries:
            while True:
                choice = input("Enter the number of the article you want to retrieve (or 'q' to quit): ")
                if choice.lower() == 'q':
                    break
                try:
                    index = int(choice) - 1
                    if 0 <= index < len(entries):
                        doi = entries[index].get('prism:doi')
                        if doi:
                            article = retrieve_article(doi)
                            if article:
                                display_article(article)
                        else:
                            print("No DOI available for this article.")
                    else:
                        print("Invalid article number. Please try again.")
                except ValueError:
                    print("Invalid input. Please enter a number or 'q' to quit.")

if __name__ == "__main__":
    main()