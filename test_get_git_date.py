# %%
import subprocess


def get_latest_commit_date() -> str:
    """現在のツリーで直近のコミット日時を文字列として取得する

    Returns
    -------
    str
        直近のコミット日時
    """
    git_command = ["git", "log", "-1", "--format='%cI'"]

    proc = subprocess.run(
        git_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    latest_commit_date_text = proc.stdout.decode("utf-8")

    return latest_commit_date_text


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
    print(get_latest_commit_date())

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
            text_about_staging = "and Not staged"
        elif proc3_text[0:2] == "MM":
            text_about_staging = "and Partially staged"
        elif proc3_text[0:2] == "M ":
            text_about_staging = "and All staged"
        else:
            text_about_staging = "and Can't read current status. Check commands for subprocess."

    print(text_about_change + text_about_staging)
