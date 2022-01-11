import subprocess

#
# command = 'job_output=$(bsub -J "{name}[1-{N}]" -n 1 -r -W 0{runtime}:05 -oo {log} {instruction})'.format(
#     name=NAME,
#     N=str(N_SIMS),
#     runtime=RUNTIME,
#     log=JI_LOG,
#     instruction=INSTRUCTION)
#
# print(command)
# system(command)

# subprocess.run(['var=3'])
output = subprocess.check_output(['echo', 'hello'])
print(output)

out = output.strip()

print(out)