"""
Get paper texts from urls.
"""
import functools
import requests
from bs4 import BeautifulSoup

PAPER_SOUP_CACHE_TIMEOUT = 10

@functools.lru_cache(maxsize=PAPER_SOUP_CACHE_TIMEOUT)
def get_soup(url):
    """
    Gets all the text from the html of a url of a paper.
    Returns the prettified text as a string.
    Requires the following libraries:
    * requests
    * BeautifulSoup from bs4
    """
    # Get the response to the request to the url
    request = requests.get(url, allow_redirects=True)
    # Get the beautiful soup
    soup = BeautifulSoup(request.content, 'html.parser')
    return soup

@functools.lru_cache(maxsize=PAPER_SOUP_CACHE_TIMEOUT)
def get_html(url: str) -> str:
    """
    Return a prettified version of the html.
    """
    soup = get_soup(url)
    return soup.prettify()

@functools.lru_cache(maxsize=PAPER_SOUP_CACHE_TIMEOUT)
def get_text(url: str) -> str:
    """
    Convert the html in the beautiful soup into text (type string)
    """
    soup = get_soup(url)
    return soup.get_text()
