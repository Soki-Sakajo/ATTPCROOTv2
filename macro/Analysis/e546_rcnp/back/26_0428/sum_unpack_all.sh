#!/bin/bash

NCORE=$(($(nproc)/3))

#76 matm: 7-112
#seq 7 51 | parallel -j "$NCORE" --load 70% --memfree 6G --joblog joblog_76matm.log \
#		     'nice -n 10 root -l -q "unpack_rcnp_76matm.C({})"'
#seq 53 112 | parallel -j "$NCORE" --load 70% --memfree 6G --joblog joblog_76matm.log \
#		     'nice -n 10 root -l -q "unpack_rcnp_76matm.C({})"'

seq 7 112 | parallel -j "$NCORE" --load 70% --memfree 6G --joblog joblog_76matm.log \
		     'nice -n 10 root -l -q "unpack_rcnp_76matm.C({})"'

#129 matm: 113-176
#seq 113 115 | parallel -j "$NCORE" --load 70%  --memfree 6G --joblog joblog_129matm.log \
#		       'nice -n 10 root -l -q "unpack_rcnp_129matm.C({})"'
#seq 117 176 | parallel -j "$NCORE" --load 70%  --memfree 6G --joblog joblog_129matm.log \
#		       'nice -n 10 root -l -q "unpack_rcnp_129matm.C({})"'

seq 113 176 | parallel -j "$NCORE" --load 70%  --memfree 6G --joblog joblog_129matm.log \
		       'nice -n 10 root -l -q "unpack_rcnp_129matm.C({})"'


