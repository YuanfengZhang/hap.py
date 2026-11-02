#!/bin/bash

# Simple performance and consistency test.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "PG evaluation test for ${HCVERSION} from ${HCDIR}"

HG19=${DIR}/../../example/chr21.fa

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

HAS_VCFEVAL=$(cat ${HCDIR}/../lib/hap.py/Haplo/version.py | grep has_vcfeval | perl -pe 's/.*([0-9]+)$/$1/')

if [[ $HAS_VCFEVAL != 0 ]]; then
	# run hap.py with vcfeval
	${PYTHON} ${HCDIR}/hap.py \
				 	-l chr21 \
				 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
				 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
				 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
				 	-r ${DIR}/../../example/chr21.fa \
				 	-o ${TMP_OUT}.vcfeval \
				 	--engine=vcfeval \
				 	-X \
				 	--force-interactive

	if [[ $? != 0 ]]; then
		echo "hap.py failed!"
		exit 1
	fi

	mkdir -p "${TMPDIR}/summary-legacy-evidence" || exit 1

	${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.vcfeval.summary.csv ${DIR}/../../example/happy/expected.vcfeval.summary.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.vcfeval.summary.csv").legacy-compare.log" 2>&1
	LEGACY_SUMMARY_STATUS=$?
	printf '%s\n' "$LEGACY_SUMMARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.vcfeval.summary.csv").legacy-compare.exit-code"
	${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.vcfeval.summary.csv ${DIR}/../../example/happy/expected.vcfeval.summary.csv --schema ga4gh-summary-v2
	if [[ $? != 0 ]]; then
		echo "vcfeval summary differs! -- diff ${TMP_OUT}.vcfeval.summary.csv ${DIR}/../../example/happy/expected.vcfeval.summary.csv"
		exit 1
	fi	
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT} \
			 	-X \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# This script checks if the summary precision / recall figures have changed significantly
mkdir -p "${TMPDIR}/summary-legacy-evidence" || exit 1
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected.summary.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.log" 2>&1
LEGACY_SUMMARY_STATUS=$?
printf '%s\n' "$LEGACY_SUMMARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.exit-code"
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected.summary.csv --schema ga4gh-summary-v2
if [[ $? != 0 ]]; then
	echo "All summary differs! -- diff ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected.summary.csv"
	exit 1
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT}.pass \
			 	-X --pass-only \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

mkdir -p "${TMPDIR}/summary-legacy-evidence" || exit 1

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/happy/expected.pass.summary.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.pass.summary.csv").legacy-compare.log" 2>&1
LEGACY_SUMMARY_STATUS=$?
printf '%s\n' "$LEGACY_SUMMARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.pass.summary.csv").legacy-compare.exit-code"
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/happy/expected.pass.summary.csv --schema ga4gh-summary-v2
if [[ $? != 0 ]]; then
	echo "PASS summary differs! -- diff ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/happy/expected.pass.summary.csv"
	exit 1
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT}.unhappy \
			 	-X --unhappy \
			 	--roc INFO.VQSLOD \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

mkdir -p "${TMPDIR}/summary-legacy-evidence" || exit 1

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.unhappy.summary.csv ${DIR}/../../example/happy/expected.unhappy.summary.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.unhappy.summary.csv").legacy-compare.log" 2>&1
LEGACY_SUMMARY_STATUS=$?
printf '%s\n' "$LEGACY_SUMMARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.unhappy.summary.csv").legacy-compare.exit-code"
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.unhappy.summary.csv ${DIR}/../../example/happy/expected.unhappy.summary.csv --schema ga4gh-summary-v2 --qq-field INFO.VQSLOD
if [[ $? != 0 ]]; then
	echo "PASS summary differs! -- diff ${TMP_OUT}.unhappy.summary.csv ${DIR}/../../example/happy/expected.unhappy.summary.csv"
	exit 1
fi

rm -f ${TMP_OUT}*
echo "PG quantification test SUCCEEDED!"

