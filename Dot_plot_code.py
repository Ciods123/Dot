import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

df = pd.read_excel(r"C:\Users\admin\Downloads\SLK_Autophagy (1).xlsx")

def scale_number(unscaled, to_min, to_max, from_min, from_max):
    return (to_max - to_min) * (unscaled - from_min) / (from_max - from_min) + to_min

def scale_list(l, to_min, to_max):
    return [scale_number(i, to_min, to_max, min(l), max(l)) for i in l]

lg2cut = 0.37
min_val = df['log2fc_log10fc'].min()
max_val = df['log2fc_log10fc'].max()

v = scale_list([min_val, (min_val + (-lg2cut)) / 2, -lg2cut, 0, lg2cut, (lg2cut + max_val) / 2, max_val], 0, 1)
print(v)

c = ["#00379D", "#295CBB", "#5386E5", "#E3F2FD", "#F96262", "#E33838", "#D10000"]
l = list(zip(v, c))
print(l)

cmap = LinearSegmentedColormap.from_list('rg', l, N=256)

exp_conditions = df['Cancer type']
proteins = df['gene_site']
log2fc_log10fc = df['log2fc_log10fc']
# bubble_sizes = df['bubble_size']  # Assuming there is a column 'bubble_size' for bubble sizes
# bubble_sizes = bubble_sizes.abs() * 100
bubble_sizes = log2fc_log10fc.abs() * 350

plt.figure(figsize=(30, 65)) 
ax = plt.gca()
ax.set_facecolor('#FFF')

# Adjust the s parameter with 'bubble_sizes' to change the bubble size
# scatter_plot = plt.scatter(exp_conditions, proteins, s=bubble_sizes, c=log2fc_log10fc, cmap=cmap, marker='o')
scatter_plot = plt.scatter(exp_conditions, proteins, s=bubble_sizes, c=log2fc_log10fc, cmap=cmap, marker='o')

colorbar = plt.colorbar(scatter_plot,aspect = 20,shrink=0.3,pad=0.13)
colorbar.ax.tick_params(labelsize=50)  
colorbar.set_label('log2fc_log10fc', fontsize=45)


expression_legend = plt.legend(handles=[
    Patch(facecolor='#D10000', label='up-regulation', edgecolor='black'),
    Patch(facecolor='#E3F2FD', label='no-expression', edgecolor='black'),
    Patch(facecolor='#00379D', label='down-regulation', edgecolor='black')
], title='Expression', bbox_to_anchor=(0.8, 0.87), fancybox=True, framealpha=1.0, bbox_transform=plt.gcf().transFigure,
loc='upper left', fontsize=40, title_fontsize=30)


bubble_size_legend_handles = []
for size in [700, 2000, 6000]:  # DOWN
    bubble_size_legend_handles.append(plt.scatter([], [], s=size, color='#00379D'))

bubble_size_legend = plt.legend(bubble_size_legend_handles, ['0 -> -5', '-5 -> -10', '-10 -> -15'], title='DOWN % Exp', bbox_to_anchor=(0.75, 0.9), fancybox=True, framealpha=1.0, bbox_transform=plt.gcf().transFigure,
loc='upper left', fontsize=50, title_fontsize=45)


bubble_size_legend_handles_1 = []

for size in [700, 1800, 6000]:  # UP
    bubble_size_legend_handles_1.append(plt.scatter([], [], s=size, color='#D10000'))

bubble_size_legend_1 = plt.legend(bubble_size_legend_handles_1, ['0 -> 3', '3 -> 6', '6 -> 9'], title='UP % Exp', bbox_to_anchor=(0.75, 0.98), fancybox=True, framealpha=1.0, bbox_transform=plt.gcf().transFigure,
loc='upper left', fontsize=50, title_fontsize=45)



# plt.gca().add_artist(expression_legend)
plt.gca().add_artist(bubble_size_legend) 
plt.gca().add_artist(bubble_size_legend_1) 




plt.xlabel('Celltype disease', fontsize=45)
plt.ylabel('Gene_site', fontsize=45)
plt.title('SLK', fontsize=50)

plt.xticks(rotation=90)  
plt.tick_params(axis='x', labelsize=35)  
plt.tick_params(axis='y', labelsize=28) 
plt.tight_layout() 

plt.savefig("SLK_autophagy_plot_new1.svg", format="svg", bbox_inches='tight')
plt.show()



