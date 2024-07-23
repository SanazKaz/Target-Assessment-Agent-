import requests
import json
from dotenv import load_dotenv
import os

load_dotenv()

API_KEY = os.getenv('SCIENCE_DIRECT_API_KEY')
BASE_URL = 'https://api.elsevier.com/content/article/doi/'

def get_article(doi):
    url = BASE_URL + doi
    params = {
        'apiKey': API_KEY
    }
    headers = {
        'Accept': 'application/json'
    }
    response = requests.get(url, params=params, headers=headers)
    limit_time = response.headers.get('X-RateLimit-Limit') # Get the  quote limit
    time_remaining = response.headers.get('X-RateLimit-Remaining') # Get the remaining quote limit
    reset_time = response.headers.get('X-RateLimit-Reset') # Get the time when the quote limit will reset
    
    print(f"Rate limit reset time: {reset_time}")
    print(f"Rate limit: {limit_time}")
    print(f"Time remaining: {time_remaining}")

    if response.status_code == 200:
        print(f"Status code: {response.status_code}")
        return response.json()
    else:
        print(f"Error: {response.status_code}")
        print(response.text)
        return None

def main():
    doi = input("Enter your DOI: ")
    results = get_article(doi)
    
    if results:
        print(json.dumps(results, indent=2))
    else:
        print("Unable to retrieve article.")

if __name__ == "__main__":
    main()