from scipy.stats import f_oneway, mannwhitneyu
import itertools


def statistical_test_same_well(ax, data, significance_level=0.05):
    pairs = list(itertools.combinations(range(len(data)), 2))
    y_max = max(max(sublist) for sublist in data)
    height = y_max * 0.05  # Adjust height of the significance marker
    for pair in pairs:
        group1 = data[pair[0]]
        group2 = data[pair[1]]
        # Perform Mann-Whitney U test
        stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
        
        if p < significance_level:
            x1, x2 = pair
            y, h = y_max + height, height
            # Determine the significance level
            if p < 0.001:
                significance = '***'
            elif p < 0.01:
                significance = '**'
            elif p < 0.05:
                significance = '*'
            else:
                significance = ''
            # Plot significance line and asterisk
            if significance:
                ax.plot([x1 + 1, x1 + 1, x2 + 1, x2 + 1], [y, y + h, y + h, y], lw=1.5, c='black')
                ax.text((x1 + x2 + 2) * .5, y + h, significance, ha='center', va='bottom', color='black')
            # Increase y_max for the next line
            y_max += height * 2