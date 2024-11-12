import pandas as pd

def csv_to_dataflame(FILE)-> None:
    df = pd.read_csv(FILE)
    all_name_list = df.iloc[0:,0].unique().tolist()
    return df, all_name_list
