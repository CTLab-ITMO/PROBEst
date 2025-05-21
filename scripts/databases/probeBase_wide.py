import pandas as pd

# Read the CSV file
df = pd.read_csv('data/databases/open/probeBase.csv', header=None)

# The first column contains the attribute names, second column contains values, 
# and third column contains the id
df.columns = ['Attribute', 'Value', 'id']

# Pivot the dataframe
df_wide = df.pivot(index='id', columns='Attribute', values='Value')

# Reset the index to make 'id' a regular column
df_wide.reset_index(inplace=True)

# Save to CSV
df_wide.to_csv('data/databases/open/probeBase_formatted.csv', index=False)