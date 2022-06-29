# %%
import subprocess
from subprocess import PIPE

if __name__ == "__main__":
    proc = subprocess.run(["git", "show", "HEAD"], stdout=PIPE, stderr=PIPE, shell=True)
    print(proc.stdout.decode("cp932"))
