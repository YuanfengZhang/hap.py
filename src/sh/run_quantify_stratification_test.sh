#!/bin/bash

# Quantification and GA4GH intermediate file format compliance
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

# Independent fix contract; does not replace the historical golden comparisons below.
CONF_UNION_ROOT=$(mktemp -d "${TMPDIR}/happy-conf-union.XXXXXXXXXX") || exit 1
"${PYTHON}" -I "${DIR}/run_conf_union_regression.py" \
    --bin-dir "${HCDIR}" --work-dir "${CONF_UNION_ROOT}/results" || exit 1

echo "Quantify stratified counting test for ${HCVERSION} from ${HCDIR}"

HG38=${DIR}/../../example/happy/hg38.chr21.fa

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py on NA12878 test file
# note this only works if we don't change the confident regions
# based on the truthset
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_hg38-chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878-GATK3-chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_hg38-chr21.bed.gz \
			 	-r ${HG38} \
			 	-o ${TMP_OUT} \
                --stratification ${DIR}/../../example/happy/stratification.tsv \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# This script checks if the summary precision / recall figures have changed significantly
mkdir -p "${TMPDIR}/summary-legacy-evidence" || exit 1
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-stratified.summary.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.log" 2>&1
LEGACY_SUMMARY_STATUS=$?
printf '%s\n' "$LEGACY_SUMMARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.exit-code"
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-stratified.summary.csv --schema ga4gh-summary-v2 > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").pre-boundary.log" 2>&1
PRE_BOUNDARY_STATUS=$?
printf '%s\n' "$PRE_BOUNDARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").pre-boundary.exit-code"
# Interval merging is fixed; compare directly with repository golden data without the six-record boundary adapter.
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-stratified.summary.csv --schema ga4gh-summary-v2
if [[ $? != 0 ]]; then
	echo "All summary differs! -- diff ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-stratified.summary.csv"
	exit 1
fi

# This script checks if the extended precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_extended.py ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/expected-stratified.extended.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.extended.csv").legacy-compare.log" 2>&1
LEGACY_EXTENDED_STATUS=$?
printf '%s\n' "$LEGACY_EXTENDED_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.extended.csv").legacy-compare.exit-code"
${PYTHON} ${DIR}/compare_extended.py ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/expected-stratified.extended.csv --schema ga4gh-metrics-v2
if [[ $? != 0 ]]; then
	echo "All extended differs! -- diff ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/expected-stratified.extended.csv"
	exit 1
fi

rm -rf ${TMP_OUT}.*
