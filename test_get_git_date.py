# %%
import subprocess
import importlib

import observation_simulation_myclass as osm


if __name__ == "__main__":
    importlib.reload(osm)

    print("--- Get all ---")
    proc1 = subprocess.run(
        ["git", "log", "-1"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)
    proc1_text = proc1.stdout.decode("utf-8")
    print(proc1_text)

    print("--- Get only date---")
    print(osm.get_latest_commit_date())

    print("--- Get change status---")
    print(osm.have_some_change_in_git_status())
