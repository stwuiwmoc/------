# %%
import subprocess

if __name__ == "__main__":

    print("--- Get all ---")
    proc1 = subprocess.run(
        ["git", "log", "-1"],
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

    print("--- Get change status---")
    proc3 = subprocess.run(
        ["git", "status"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    proc3_text = proc3.stdout.decode("utf-8")
    print(proc3_text)

    if len(proc3_text) == 0:
        print("Not changed")
    elif proc3_text[1] == "M":
        print("Not staged")
    else:
        print("Staged")
