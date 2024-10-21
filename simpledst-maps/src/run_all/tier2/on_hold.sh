#!/bin/bash

# Get the node name (which is the job name)
job_name=$1

#Get the job id 
job_id=$(condor_q -format "%d\n" ClusterId -constraint "DAGNodeName == \"$job_name\"")

# Check if job exists
if ! condor_q $job_id &>/dev/null; then
    echo "Job $job_id does not exist or is no longer valid." 
    exit 1
fi

# Get the actual memory used by the job before it was held
hold_reason=$(condor_q -long $job_id | grep "HoldReason")
actual_memory_used=$(echo $hold_reason | grep -oP 'used \K[0-9]+(?=MB)')

# Check the current memory request
current_memory=$(condor_q -long $job_id | grep "RequestMemory" | awk '{print $3}')

# Check if memory retrieval was successful
if [ -z "$current_memory" ] || [ -z "$actual_memory_used" ]; then
    echo "Could not retrieve current memory for job $job_id."
    exit 1
fi

# Increase memory to actual amount + 1 GB 
new_memory=$(($actual_memory_used + 1000))

# Set a maximum memory limit
max_memory=10000
if [ "$new_memory" -gt "$max_memory" ]; then
    echo "Memory request for job $job_id exceeds limit ($max_memory MB). Aborting job." 
    condor_rm $job_id
    exit 1
fi

# Update the memory request for the held job
if ! condor_qedit $job_id RequestMemory $new_memory; then
    echo "Failed to update memory for job $job_id."
    exit 1
fi

# Release the held job
if ! condor_release $job_id; then
    echo "Failed to release job $job_id."
    exit 1
fi

echo "Job $job_id memory updated to $new_memory MB and released." 
