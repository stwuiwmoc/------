# %%
import subprocess

if __name__ == "__main__":
    print("--- start ---")
    proc1 = subprocess.run(
        ["git", "show", "HEAD"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)
    print(proc1.stdout.decode("utf-8"))

    print("--- start ---")
    proc2 = subprocess.run(
        ["git", "log"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)
    print(proc2.stdout.decode("utf-8"))
