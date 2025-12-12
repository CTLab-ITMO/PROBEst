import requests
from bs4 import BeautifulSoup
import pandas as pd


def parse_probebase_page(url: str) -> pd.Series:
    """Parse one page from ProbeBase database to the uniform format

    Parameters
    ----------
    url : str
        URL string, path to the probebase page 

    Returns
    -------
    table : pd.DataFrame
        parsed table from probebase
    """

    # Download html and get table

    response = requests.get(url)

    # Response checking
    if response.status_code == 200:
        # Parse the HTML content using BeautifulSoup
        soup = BeautifulSoup(response.content, 'html.parser')

        # if table is empty
        if soup.find_all('tr') == []:  # CORRECT!!!!!!!!!!!!!!!!!!!
            return Warning('Page without probe data.frame or parsing problems')

        # Create a list to hold the rows of the table
        table_data = []

        # Loop through each row in the table
        for row in soup.find_all('tr'):
            # print(row, "###############")
            cells = row.find_all(['td', 'th', 'value'])
            row_data = [cell.get_text(strip=True) for cell in cells]
            table_data.append(row_data)

        # Convert the list of rows into a DataFrame for easier manipulation
        df = pd.DataFrame(table_data)
        
        # Check if we have enough columns
        if df.shape[1] < 2:
            return Warning('Table has insufficient columns')
            
        df.iloc[0,0] = 'Test'
        df.iloc[1,0] = 'Name'

        df2 = df.iloc[:,1]
        df2.index = df.iloc[:,0]

        # Display the result
        return df2

    else:
        return ImportWarning(f"Failed to retrieve the page. Status code: {response.status_code}")