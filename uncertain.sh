#!/bin/bash

# \todo: keep CONAN_VERSION updated, test thoroughly whenever you do, leave this "todo" here
CONAN_VERSION=1.59
GCOVR_VERSION=3.4
COVERAGE_INPUT_DIR=build
COVERAGE_OUTPUT_DIR=coverage_reports
STATIC_CHECKS_DIR=static-checks

# Fail if any command fails
set -e
set -o pipefail

# Print each command as it executes
QUIET="--quiet"
if [ -n "$VERBOSE" ]; then
  set -o xtrace
  QUIET=""
fi

# This script should work no matter where you call it from.
cd "$(dirname "$0")"

FAILURE=1
SUCCESS=0

# Check if Darwin or Linux
case "$OSTYPE" in
  darwin*)  echo "Running Darwin: ok" ;;
  linux*)   echo "Running Linux: ok" ;;
  *)
    echo "Error: This script only supports 'Darwin' or 'Linux'. You have $OSTYPE."
    exit $FAILURE
    ;;
esac

####################
# HELPER FUNCTIONS #
####################

print_help() {
    cat <<EOF

The following options are available:

EOF
}

ensure_not_root() {
  if [ "$EUID" -eq 0 ] && [ -z "$FORCED_ROOT" ]; then
    echo "Please do not run with root privileges!"
    exit $FAILURE
  fi
}

get_num_cpus() {
  # build with 1 less than total number of CPUS, minimum 1
  NUM_CPUS=$(cat /proc/cpuinfo | grep -c processor)
  let NUM_CPUS-=1
  if [ "$NUM_CPUS" -lt "1" ]; then
    NUM_CPUS=1
  fi
  echo ${NUM_CPUS}
}

clean_dir() {
  dir_name=$1
  if [ -d "$dir_name" ]; then
    echo "Removing $dir_name"
    rm -rf "$dir_name"
    return $SUCCESS
  elif [ -f "$dir_name" ]; then
    echo "File with this name already exists, not a directory."
    return $FAILURE
  fi
}

create_clean_directory() {
  dir_name=$1
  clean_dir "$dir_name"
  if mkdir -p "$dir_name"; then
    echo "Clean directory created: $dir_name"
    return $SUCCESS
  else
    echo "Creating directory failed: $dir_name"
    return $FAILURE
  fi
}

install_linux() {
  sudo apt-get update
  sudo apt-get install -y cmake python3-pip python3-setuptools clang-tidy cppcheck
  sudo pip3 install -U pip
}

install_pip() {
  pip3 install install conan==$CONAN_VERSION
  pip3 install install gcovr==$GCOVR_VERSION
  pip3 install install codecov cpp-coveralls clang-html
  source "${HOME}/.profile"
  conan profile new --detect default
  conan profile update settings.compiler.libcxx=libstdc++11 default
}

run_cppcheck() {
  CPPCHECK_OUTPUT_DIR=${STATIC_CHECKS_DIR}/cppcheck
  if [ "$1" == "html" ]; then
    create_clean_directory ${CPPCHECK_OUTPUT_DIR}
    xml_opts="--xml --output-file=${CPPCHECK_OUTPUT_DIR}/report.xml"
  fi

  cppcheck --language=c++ --enable=all --inconclusive --force --inline-suppr ${QUIET} \
           --project=build/compile_commands.json \
           --suppress=missingIncludeSystem \
           --suppress=*:*/gtest/* \
           --addon=misc \
           ${xml_opts} \
           .

#           --addon=misra.py \

  if [ "$1" == "html" ]; then
    cppcheck-htmlreport --file=${CPPCHECK_OUTPUT_DIR}/report.xml \
                        --title="Uncertain lib" \
                        --report-dir=${CPPCHECK_OUTPUT_DIR} --source-dir=.
    python -m webbrowser "${CPPCHECK_OUTPUT_DIR}/index.html"
  fi
}

run_clang_tidy() {
  j_opt=$2

  CLANG_TIDY_EXEC="run-clang-tidy"
  CLANG_TIDY_VERSION=$(echo "$(clang-tidy --version | sed -n 1p)" | awk -F[" ".] '{print $4}')

  if [ "$1" == "html" ]; then
    TIDY_OUTPUT_DIR=${STATIC_CHECKS_DIR}/clang-tidy
    TIDY_FILE=${TIDY_OUTPUT_DIR}/logfile.log
    create_clean_directory ${TIDY_OUTPUT_DIR}
  else
    TIDY_FILE="/dev/stdout"
  fi

  files=$(find . -regextype posix-extended -regex '.*\.(cpp|hpp|h|tpp)' -not -path "*build*")

  echo "Will run $CLANG_TIDY_EXEC $j_opt with version=${CLANG_TIDY_VERSION}"
  echo -e "Will analyze the following files:\n${files}"

  $CLANG_TIDY_EXEC -quiet ${j_opt} -p build ${files} \
   -header-filter='^.*\/(source|examples|tests)\/.*\.(hpp|cpp|h|tpp)$'\
   > ${TIDY_FILE}

  if [ "$1" == "html" ]; then
    sed -i "s@$(pwd)@@g" ${TIDY_FILE}
    clang-tidy-html ${TIDY_FILE} -o ${TIDY_OUTPUT_DIR}/clang.html \
    -d "https://releases.llvm.org/${CLANG_TIDY_VERSION}.0.0/tools/clang/tools/extra/docs/clang-tidy/checks/list.html"
    python -m webbrowser "${TIDY_OUTPUT_DIR}/clang.html"
  fi
}

generate_coverage_reports() {
  echo "Generating test coverage reports..."

  create_clean_directory ${COVERAGE_OUTPUT_DIR}

  lcov ${QUIET} --directory "$COVERAGE_INPUT_DIR" --capture \
       --output-file "$COVERAGE_OUTPUT_DIR/coverage.info"

  lcov ${QUIET} --remove "$COVERAGE_OUTPUT_DIR/coverage.info" \
       --output-file "$COVERAGE_OUTPUT_DIR/coverage_trimmed.info" \
       "*/tests/*" \
       "*gtest*" \
       "*spdlog*" \
       "*fmt*" \
       "/usr/include*"

  rm "$COVERAGE_OUTPUT_DIR/coverage.info"
  mv "$COVERAGE_OUTPUT_DIR/coverage_trimmed.info" "$COVERAGE_OUTPUT_DIR/coverage.info"
}

view_coverage() {
  genhtml ${QUIET} "$COVERAGE_OUTPUT_DIR/coverage.info" \
      --output-directory "$COVERAGE_OUTPUT_DIR"

  echo "Coverage reports generated at '$COVERAGE_OUTPUT_DIR/index.html'"
  echo "   You may open it in browser with 'python3 -m webbrowser ${COVERAGE_OUTPUT_DIR}/index.html'"

  python -m webbrowser "${COVERAGE_OUTPUT_DIR}/index.html"
}

upload_coverage_reports() {
  echo "Uploading coverage reports to Codecov"

  curl -Os https://uploader.codecov.io/latest/linux/codecov
  chmod +x codecov
  ./codecov
  rm codecov
}

########
# HELP #
########

if [ "$1" == "help" ] || [ "$1" == "-h" ]; then
  print_help
  exit $SUCCESS

###########
# INSTALL #
###########
elif [ "$1" == "install" ]; then
  case "$OSTYPE" in
    linux*)
      ensure_not_root
      install_linux
      install_pip
      exit $SUCCESS
      ;;
    *)
      echo "Error: No installation scripts for $OSTYPE."
      exit $FAILURE
      ;;
  esac

#########
# CLEAN #
#########
elif [ "$1" == "clean" ]; then
  clean_dir build
  clean_dir "$COVERAGE_OUTPUT_DIR"
  exit $SUCCESS

#########
# BUILD #
#########
elif [ "$1" == "build" ]; then

  config_type="Release"
  if [ "$2" == "--debug" ] || [ "$3" == "--debug" ] || [ "$4" == "--debug" ]; then
    config_type="Debug"
  fi

  j_opt=""
  if [ "$2" == "-j" ] || [ "$3" == "-j" ] || [ "$4" == "-j" ]; then
    j_opt="-j$(eval get_num_cpus)"
  fi

  checks_opt="-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
  if [ "$2" == "--no-checks" ] || [ "$3" == "--no-checks" ] || [ "$4" == "--no-checks" ]; then
    checks_opt=""
  fi

  ensure_not_root
  create_clean_directory build
  pushd build
  cmake -DCMAKE_BUILD_TYPE=${config_type} ${checks_opt} ..
  make everything ${j_opt}
  popd

  exit $SUCCESS

#########
# CHECK #
#########
elif [ "$1" == "check" ]; then
  if [ "$2" == "-j" ] || [ "$3" == "-j" ]; then
    j_opt="-j$(eval get_num_cpus)"
  fi

  if [ "$2" == "--html" ] || [ "$3" == "--html" ]; then
    run_clang_tidy html $j_opt
    run_cppcheck html
  else
    run_clang_tidy nohtml $j_opt
    run_cppcheck
  fi



########
# TEST #
########
elif [ "$1" == "test" ]; then

  j_opt=""
  if [ "$2" == "-j" ]; then
    j_opt="-j$(eval get_num_cpus)"
  fi

  ensure_not_root
  create_clean_directory build
  pushd build
  cmake -DCMAKE_BUILD_TYPE=Debug -DCOV=1 ..
  make run_tests ${j_opt}
  make examples ${j_opt}
  ./bin/ms_demo > null
  ./bin/full_demo > null
  popd

  exit $SUCCESS

############
# Coverage #
############

elif [ "$1" == "cov" ]; then
  ensure_not_root
  generate_coverage_reports
  if [ "$2" == "--upload" ]; then
    upload_coverage_reports
  else
    view_coverage
  fi
  exit $SUCCESS

################
# ERROR & HELP #
################
else
  echo "No valid options provided :\("
  print_help
  exit $FAILURE
fi
