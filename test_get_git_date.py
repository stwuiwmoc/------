# %%
import subprocess

if __name__ == "__main__":
    proc1 = subprocess.run(
        ["git", "show", "HEAD"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)
    proc1_text = proc1.stdout.decode("utf-8")

    target = "Date:   "
    idx = proc1_text.find(target)
    extracted_text = proc1_text[idx + 12:idx + 38]
    print(extracted_text)
