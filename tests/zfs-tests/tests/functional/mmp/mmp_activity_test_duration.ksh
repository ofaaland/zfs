#!/bin/ksh -p
#
# CDDL HEADER START
#
# This file and its contents are supplied under the terms of the
# Common Development and Distribution License ("CDDL"), version 1.0.
# You may only use this file in accordance with the terms of version
# 1.0 of the CDDL.
#
# A full copy of the text of the CDDL should have accompanied this
# source.  A copy of the CDDL is also available via the Internet at
# http://www.illumos.org/license/CDDL.
#
# CDDL HEADER END
#

#
# Copyright (c) 2017 by Lawrence Livermore National Security, LLC.
#

# DESCRIPTION:
#	mmp_activity_test kstat should report activity test remaining time
#
# STRATEGY:
#	1. Create a pool and enable multihost
#	2. Force export the pool and change the HOSTID
#	3. Import the pool, triggering the activity test
#	4. Verify /proc/spl/kstat/zfs/<poolname>/mmp_activity_test contents
#

. $STF_SUITE/include/libtest.shlib
. $STF_SUITE/tests/functional/mmp/mmp.cfg
. $STF_SUITE/tests/functional/mmp/mmp.kshlib

verify_runnable "both"

function cleanup
{
	# Let the import finish, then clean up
	sleep 15
	default_cleanup_noexit
	log_must mmp_clear_hostid
}

ACT_TEST_LEFT_KSTAT="/proc/spl/kstat/zfs/\$import/activity_test"

log_assert "mmp_activity_test kstat reports activity test remaining time"
log_onexit cleanup

# 1. Create a zpool and enable multihost
log_must mmp_set_hostid $HOSTID1
default_setup_noexit $DISK
log_must zpool set multihost=on $TESTPOOL

# 2. Force export the pool and change the hostid
log_must zpool export -F $TESTPOOL
log_must mmp_clear_hostid
log_must mmp_set_hostid $HOSTID2

# 3. Import the pool in the background
zpool import -f $TESTPOOL >/dev/null 2>&1 &

# 4. Wait for the kernel to start the import, then verify the remaining time
# Expected time for the import test is 10-12.5 seconds
# Expected time for import process is < 1 second
# So after waiting 2 seconds to get into the activity test,
#   8s <= remaining <= 11.5s
EXPECTED_MIN=8000
EXPECTED_MAX=11500

sleep 2
log_must pgrep zpool
log_must cat $ACT_TEST_LEFT_KSTAT

remaining_ms=$(($(awk '{print $1}' $ACT_TEST_LEFT_KSTAT) / 1000000))
if [ $remaining_ms -lt $EXPECTED_MIN ] || [ $remaining_ms -gt $EXPECTED_MAX ]; then
	log_fail "ERROR: reported $remaining_ms ms not within [$EXPECTED_MIN-$EXPECTED_MAX]"
fi

log_pass "mmp_activity_test kstat reports activity test remaining time"
