
#Import phases
#Search
#Tryimport
#ActivityCheck
#Import
#MMPThreadStarts
#
#If MachineB gets to ActivityCheck after MachineA has written labels out, a full activity check will occur because pool state is ACTIVE.
#
#Problem occurs when MachineB gets to ActivityCheck while MachineB is also there, or while it is in Import.
#
#So if ActivityCheck on either machine is longer than time from ActivityCheck to MMPThreadStarts on the other machine, the dual imports will be detected.
#
#How long is that period of time?

cat /proc/spl/kstat/zfs/dbgmsg | awk '
	/pool .[$]import. import entered spa_activity.* at / {tryimport_check=$NF}
	/pool .[^$].* import entered spa_activity.* at / {import_check=$NF}
	/pool .* import entered mmp_thread.* at / {mmp_thread=$NF}
	END { print "tryimport_check",tryimport_check,"import_check",import_check,"mmp_thread",mmp_thread}
	END { print "activity_check_interval",(import_check-tryimport_check)/1000000,"ms"}
	END { print "to_mmp_thread_interval",(mmp_thread-import_check)/1000000,"ms"}
'

# on diane vm, 3 file vdevs in /tmp
# tryimport_check 44699203654975 import_check 44699238876131 mmp_thread 44699255994801
# activity_check_interval 35.2212 ms
# to_mmp_thread_interval 17.1187 ms

# on diane vm, 3 file vdevs in /tmp referenced via loopback devices
# tryimport_check 44893629095002 import_check 44893662400853 mmp_thread 44893677606180
# activity_check_interval 33.3059 ms
# to_mmp_thread_interval 15.2053 ms

# on slag node, 70 rotating disk vdevs in SAS enclosure
# tryimport_check 3127515608005206 import_check 3127516900583564 mmp_thread 3127517514268109
# activity_check_interval 1292.58 ms
# to_mmp_thread_interval 613.685 ms

