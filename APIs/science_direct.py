import requests
import json
from dotenv import load_dotenv
import os

load_dotenv()


API_KEY = os.getenv('SCIENCE_DIRECT_API_KEY')
BASE_URL = 'https://api.elsevier.com/content/search/sciencedirect'

def search_articles(query):
    params = {
        'query': query,
        'apiKey': API_KEY
    }
    headers = {'Accept': 'application/json'}
    
    response = requests.get(BASE_URL, params=params, headers=headers)
    
    if response.status_code == 200:
        print(f"status code {response.status_code}")
        return response.json()
    else:
        print(f"Error: {response.status_code}")
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
    else:
        print("No results found or unexpected response format.")

def main():
    search_query = input("Enter your search query: ")
    results = search_articles(search_query)
    
    if results:
        display_search_results(results)

if __name__ == "__main__":
    main()