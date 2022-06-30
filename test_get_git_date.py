# %%
import subprocess

if __name__ == "__main__":

    print("--- Get all ---")
    proc1 = subprocess.run(
        ["git", "show", "HEAD"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)
    proc1_text = proc1.stdout.decode("utf-8")
    print(proc1_text)

    print("--- Get only date---")

    proc2 = subprocess.run(
        ["git", "log", "-1", "--format='%cI'"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)
    proc2_text = proc2.stdout.decode("utf-8")
    print(proc2_text)
