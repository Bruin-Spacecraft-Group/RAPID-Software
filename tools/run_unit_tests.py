import os
import shutil
import argparse
import subprocess
import sys

import deploy_to_usb


def RED(text):
    return "\033[91m" + text + "\033[0m"


def GREEN(text):
    return "\033[92m" + text + "\033[0m"


def YELLOW(text):
    return "\033[93m" + text + "\033[0m"


def CYAN(text):
    return "\033[96m" + text + "\033[0m"


exit_code_names = {
    0: "OK",
    1: "TESTS_FAILED",
    2: "INTERRUPTED",
    3: "INTERNAL_ERROR",
    4: "USAGE_ERROR",
    5: "NO_TESTS_COLLECTED",
}


parser = argparse.ArgumentParser(
    description="Run this script to deploy each artifact in the artifacts folder "
    + "to a temporary folder and run the unit tests which apply to that artifact."
)
args = parser.parse_args()

artifacts_dir = os.path.join(".", "artifacts")
tests_passed = True

for app in os.listdir(artifacts_dir):
    if app.startswith(".DS_Store"):
        print(f"Not generating unit tests for artifact {app}: .DS_Store directory.")
        continue
    if app.endswith("_testapp"):
        print(f"Not generating unit tests for artifact {app}: already a test app")
        continue

    print(CYAN(f"Generating tests for artifact {app}..."))
    test_app_dir = deploy_to_usb.deploy_with_settings(app, None, True, True)

    print(CYAN(f"Running pytest for deployment of {app}..."))
    test_env = os.environ.copy()
    test_env["RAPID-0-Software_test-config"] = os.path.abspath(
        os.path.join(".", "config")
    )
    result = subprocess.run(
        [sys.executable, "-P", "-m", "pytest"], cwd=test_app_dir, env=test_env
    ).returncode
    if result != 0:
        tests_passed = False
        print(
            f"Pytest Exit Code {result}: {exit_code_names.get(result, "UNKNOWN_EXIT_CODE")}"
        )

exit(0 if tests_passed else 1)
