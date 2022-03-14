import os
os.system('export NCBI_API_KEY="447726e72ee80490c9d97644aa34d3e27508"')  # set metapub key
from papers2ids import Papers
papers = Papers(selected_date="2022_03_08")
papers()
