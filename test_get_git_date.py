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
        ["git", "status", "--short"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    proc3_text = proc3.stdout.decode("utf-8")
    print(proc3_text)

    if len(proc3_text) == 0:
        text_about_change = "Not changed"
        text_about_staging = ""

    else:
        text_about_change = "Changed"

        if proc3_text[0:2] == " M":
            text_about_staging = "Not staged"
        elif proc3_text[0:2] == "MM":
            text_about_staging = "Partially staged"
        elif proc3_text[0:2] == "M ":
            text_about_staging = "All staged"
        else:
            text_about_staging = "Can't read current status. Check commands for subprocess."

    print(text_about_change + " and " + text_about_staging)
