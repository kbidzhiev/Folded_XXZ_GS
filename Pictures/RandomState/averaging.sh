
#!/bin/bash
 


paste realization_1/Sz_profile.dat realization_2/Sz_profile.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp1.dat
paste realization_3/Sz_profile.dat realization_4/Sz_profile.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp2.dat
paste realization_5/Sz_profile.dat realization_6/Sz_profile.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp3.dat
paste realization_7/Sz_profile.dat realization_8/Sz_profile.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp4.dat
paste realization_9/Sz_profile.dat realization_10/Sz_profile.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp5.dat



paste Sz_profile_tmp1.dat Sz_profile_tmp2.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp12.dat
paste Sz_profile_tmp3.dat Sz_profile_tmp4.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp34.dat

paste Sz_profile_tmp12.dat Sz_profile_tmp34.dat | awk '{ print $1 "\t"  ($2 + $6) "\t" $3 "\t" $4 ; }' > Sz_profile_tmp.dat

paste Sz_profile_tmp.dat Sz_profile_tmp5.dat | awk '{ print $1 "\t"  ($2 + $6)/10. "\t" $3 "\t" $4 ; }' > Sz_profile.dat




#   sed -i 's/search_string/replace_string/' filename

#removes 0 0 0 between blocks
sed -i 's/0\t0\t0\t0/\t/' Sz_profile.dat

#removes remain after "t=*" in each block
sed -i 's/"t=*.*"\t0/\t/' Sz_profile.dat


rm Sz_profile_tmp1.dat
rm Sz_profile_tmp2.dat
rm Sz_profile_tmp3.dat
rm Sz_profile_tmp4.dat
rm Sz_profile_tmp12.dat
rm Sz_profile_tmp34.dat
rm Sz_profile_tmp5.dat
rm Sz_profile_tmp.dat

#gnuplot Gif.gp
#eog Sz_profile.gif

#cd ~/Programs/Folded_XXZ_GS/Pictures/RandomState/realization_6/
#gnuplot  Gif.gp


