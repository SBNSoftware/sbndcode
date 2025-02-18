parentfile="parents.list"

for parent in $(cat $parentfile); do
    # if line does not start with "reco2", skip
    if [[ $parent != reco2* ]]; then
        continue
    fi
    filedir=$(samweb -e sbnd locate-file $parent)
    filedir=${filedir#enstore:}
    fullfilepath=$filedir/$parent
    echo $fullfilepath >> /exp/sbnd/data/users/munjung/job.log
    lar -c filtereventid.fcl -s $fullfilepath

    # count number of lines in filter_evts.txt
    mv filter_evts.txt filter_evts_$parent.txt
    num_lines=$(wc -l filter_evts_$parent.txt | awk '{print $1}')
    echo "number of lines in filter_evts_$parent.txt" >> /exp/sbnd/data/users/munjung/job.log
    echo $num_lines >> /exp/sbnd/data/users/munjung/job.log
    if [[ $num_lines -eq 1 ]]; then
        mv *.root empty/
        rm filter_evts_$parent.txt
    fi

    mv *.root selected/

    rm memory.db
    rm messages.log
    rm errors.log
    rm cputime.db
done

mv *list selected/.
mv *txt selected/.
mv *csv selected/.