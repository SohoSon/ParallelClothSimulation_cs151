import pandas as pd
import matplotlib.pyplot as plt

# read data
#data = pd.read_csv('spring_simulation_output.txt', names=['Threads', 'FPS', 'Total Time'], skiprows=1, dtype={'Threads': int, 'FPS': float, 'Total Time': float})
#node_based_simulation_output
data = pd.read_csv('node_based_simulation_output.txt', names=['Threads', 'FPS', 'Total Time'], skiprows=1, dtype={'Threads': int, 'FPS': float, 'Total Time': float})
# check data type
print(data.dtypes)

# calculate average FPS and total time for each thread number
average_data = data.groupby('Threads').mean().reset_index()

# plot FPS vs thread number
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(average_data['Threads'], average_data['FPS'], marker='o', linestyle='-')
plt.title('Threads vs Average FPS (Pull Model)')
plt.xlabel('Threads')
plt.ylabel('Average FPS')
plt.grid(True)

# plot total time vs thread number
plt.subplot(1, 2, 2)
plt.plot(average_data['Threads'], average_data['Total Time'], marker='o', linestyle='-', color='r')
plt.title('Threads vs Average Total Time (Pull Model)')
plt.xlabel('Threads')
plt.ylabel('Average Total Time (ms)')
plt.grid(True)

# show chart
plt.tight_layout()
plt.show()