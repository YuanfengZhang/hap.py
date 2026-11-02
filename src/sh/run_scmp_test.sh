#!/bin/bash

# Simple performance and consistency test.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

# The neighborhood-distance contract is independent of the historical golden data retained below.
DISTANCE_CONTRACT_ROOT=$(mktemp -d "${TMPDIR}/happy-distance-contract.XXXXXXXXXX") || exit 1
"${PYTHON}" -I "${DIR}/run_scmp_distance_regression.py" \
    --bin-dir "${HCDIR}" --work-dir "${DISTANCE_CONTRACT_ROOT}/results" || exit 1

echo "SCMP test for ${HCVERSION} from ${HCDIR}"


HG19=${DIR}/../../example/chr21.fa

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT} \
                --no-adjust-conf-regions \
                --no-leftshift --no-decompose \
                --engine=scmp-distance \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# This script checks if the extended precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_extended.py ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/test-distancebased.extended.csv
if [[ $? != 0 ]]; then
	echo "All extended differs! -- diff ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/test-distancebased.extended.csv"
	exit 1
fi

${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT} \
                --no-adjust-conf-regions \
                --no-leftshift --no-decompose \
                --engine=scmp-somatic \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# This script checks if the extended precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_extended.py ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/test-allelebased.extended.csv
if [[ $? != 0 ]]; then
	echo "All extended differs! -- diff ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/test-allelebased.extended.csv"
	exit 1
fi

rm -f ${TMP_OUT}*
echo "scmp test SUCCEEDED!"
