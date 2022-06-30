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


def get_current_git_status() -> str:
    """ファイルの変更箇所の有無とステージングの有無を文字列として取得する

    Returns
    -------
    str
        変更箇所の有無とステージングの有無についての文字列
    """
    git_command = ["git", "status", "--short"]

    proc = subprocess.run(
        git_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    proc_text = proc.stdout.decode("utf-8")

    if len(proc_text) == 0:
        # git status --short で出力される文字列が無い
        # つまり「ファイルの変更箇所なし」を意味する
        text_about_change = "Not changed"
        text_about_staging = ""

    else:
        text_about_change = "Changed"
        # git status --short で出力された文字列がある
        # つまり「ファイルに何かしらの変更箇所あり」を意味する

        if proc_text[0:2] == " M":
            # 何もステージングされていない
            text_about_staging = " and Not staged"

        elif proc_text[0:2] == "MM":
            # 変更箇所の一部分のみステージングされている
            text_about_staging = " and Partially staged"

        elif proc_text[0:2] == "M ":
            # 変更箇所は全てステージングされている
            text_about_staging = " and All staged"

        else:
            # 上記に該当しないので、git status で想定外の出力が出ている
            text_about_staging = " and Can't read current status. Check commands for subprocess."

    text_about_current_git_status = text_about_change + text_about_staging

    return text_about_current_git_status


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
    print(get_current_git_status())
