#!/usr/bin/env python
# Script by Allen B. Davis, shared with permission
import matplotlib.pyplot as plt
import numpy as np
import time


remaining="../../progress/all_accessions_remaining.txt"
completed="../../progress/list_complete.txt"
attempted="/../../progress/list_attempted.txt"
failed="/../../progress/list_failed.txt"


#Determine number of unique lines of each category
with open(remaining, 'r') as file:
    lines_remaining = file.readlines()

unique_lines_remaining = set(lines_remaining)
num_remaining = len(unique_lines_remaining)

with open(completed, 'r') as file:
    lines_completed = file.readlines()

# Count the number of unique lines
unique_lines_completed = set(lines_completed)
num_completed = len(unique_lines_completed)

with open(attempted, 'r') as file:
    lines_attempted = file.readlines()

# Count the number of unique lines
unique_lines_attempted = set(lines_attempted)
num_attempted = len(unique_lines_attempted)

with open(failed, 'r') as file:
    lines_failed = file.readlines()

# Count the number of unique lines
unique_lines_failed = set(lines_failed)
num_failed = len(unique_lines_failed)

# Create figure
fig = plt.figure(figsize=(6,2))

ax1 = fig.add_subplot(111) # primary ax
ax2 = ax1.twiny() # secondary ax (for the percentages)

## Values for testing/debugging
# num_completed = 523
# num_attempted = 204
# num_remaining = 140

total = num_completed + num_attempted + num_remaining + num_failed

dummy_y = 1 # dummy value for the y position of the bar plot

# Create bars. "left" indicates the starting x-value where that bar will begin (and then extend right)
b1 = ax1.barh(dummy_y, num_completed, color='forestgreen')
b2 = ax1.barh(dummy_y, num_attempted, left=num_completed, color="gold")
b3 = ax1.barh(dummy_y, num_failed, left=num_completed+num_attempted, color="darkorange")
b4 = ax1.barh(dummy_y, num_remaining, left=(num_completed+num_attempted+num_failed), color="firebrick")

# Create legend.
ax1.legend([b1, b2, b3, b4], ["Completed: %d (%.0f%%)"%(num_completed, np.round(100*num_completed/total)),
                          "In progress: %d (%.0f%%)"%(num_attempted, np.round(100*num_attempted/total)),
                          "Failed: %d (%.0f%%)"%(num_failed, np.round(100*num_failed/total)),
                          "Remaining: %d (%.0f%%)"%(num_remaining, np.round(100*num_remaining/total))],
           title="Total jobs: %d"%total, loc="right", bbox_to_anchor=(1.5,0.5))

# Hide the y axis and set x ticks
ax1.set_yticks([])
ax1.set_xlim(0, total)


# Create ticks for the percentage (upper) x axis
new_tick_locations = np.arange(0,110,10)
ax2.set_xlim(0,100)
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(["%.0f%%"%ntl for ntl in new_tick_locations])

# Label the time of creation
t = time.localtime()
t_string = time.strftime("%m/%d/%y %H:%M:%S", t)
ax1.text(s=t_string, x=total*1.12, y=dummy_y*1.35, size=12)

# Save figure in this directory
plt.savefig('../../progress/progress_bar.pdf', dpi=300, bbox_inches='tight')
