#!/usr/bin/env bash
#
# This script follows the standard template implementation from lar_ci.
# In fact, it follows it so closely, it actually is it.
#

function FATAL() {
  local -i Code="$1"
  shift
  echo "${0} -- ERROR (${Code}): $*" >&2
  exit $Code
} # FATAL()

[[ -n "$SETUP_LAR_CI" ]] || FATAL 1 "lar_ci must be set up."

TemplateScript="${LAR_CI_DIR}/bin/ci_regression_test_template.sh"
[[ -r "$TemplateScript" ]] || FATAL 2 "the template regression test script '${TemplateScript}' can't be found!"

# transfer the control to the template script
echo "Running template regression test script: '${TemplateScript}'..."
source "$TemplateScript" "$@"
