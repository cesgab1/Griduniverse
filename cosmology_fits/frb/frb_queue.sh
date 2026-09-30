cd /tmp/claude-0/-home-claude-griduniverse/631bebbd-9799-5238-817c-c2760930ea03/scratchpad
( python3 frb_h0.py fd 0.93 > frb_a.txt 2>&1; python3 frb_h0.py law > frb_b.txt 2>&1; HALO=30 python3 frb_h0.py fd 0.84 > frb_c.txt 2>&1 ) &
( HALO=80 python3 frb_h0.py fd 0.84 > frb_d.txt 2>&1; FEED=0.2 python3 frb_h0.py fd 0.84 > frb_e.txt 2>&1; FEED=0.45 python3 frb_h0.py fd 0.84 > frb_f.txt 2>&1 ) &
wait
