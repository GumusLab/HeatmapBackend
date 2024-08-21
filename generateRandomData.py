import numpy as np
import pandas as pd

# Number of rows and columns
num_rows = 20000
num_cols = 200

# Create row and column labels
row_labels = [f'row_{i+1}' for i in range(num_rows)]
col_labels = [f'col_{j+1}' for j in range(num_cols)]

# Generate random values between 1 and 5
data = np.random.randint(1, 6, size=(num_rows, num_cols))

# Create a DataFrame with the generated data and labels
df = pd.DataFrame(data, index=row_labels, columns=col_labels)

# Save the DataFrame to a CSV file
csv_file_path = 'random_dataset.csv'
df.to_csv(csv_file_path)
