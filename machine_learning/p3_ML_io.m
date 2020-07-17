%get confining pressures array
cp0=data0(:,1);
cp5=data5(:,1);
cp10=data10(:,1);
cp15=data15(:,1);
cp20=data20(:,1);
cp25=data25(:,1);
cp30=data30(:,1);
cp35=data35(:,1);
cp40=data40(:,1);
cp45=data45(:,1);
cp50=data50(:,1);
%%
mc_var0=data0(:,21);
mc_var5=data5(:,21);
mc_var15=data15(:,21);
mc_var20=data20(:,21);
mc_var25=data25(:,21);
mc_var30=data30(:,21);
mc_var40=data40(:,21);
mc_var45=data45(:,21);
mc_var50=data50(:,21);
%%
%input arrays
%1=confining pressure
%2=mc rate
%3=shear fraction
%4=energy variance
%5=moment range
%6=b value
%7=d value

i0=[cp0, mc_rate0, mc_var0, shear_fraction0, e_var0, e_rate0, moment_diff0, b_value0, d_value0];
i5=[cp5, mc_rate5, mc_var5, shear_fraction5, e_var5, e_rate5, moment_diff5, b_value5, d_value5];
i10=[cp10, mc_rate10, mc_var10, shear_fraction10, e_var10, e_rate10, moment_diff10, b_value10, d_value10];
i15=[cp15, mc_rate15, mc_var15, shear_fraction15, e_var15, e_rate15, moment_diff15, b_value15, d_value15];
i20=[cp20, mc_rate20, mc_var20, shear_fraction20, e_var20, e_rate20, moment_diff20, b_value20, d_value20];
i30=[cp30, mc_rate30, mc_var30, shear_fraction30, e_var30, e_rate30, moment_diff30, b_value30, d_value30];
i35=[cp35, mc_rate35, mc_var35, shear_fraction35, e_var35, e_rate35, moment_diff35, b_value35, d_value35];
i40=[cp40, mc_rate40, mc_var40, shear_fraction40, e_var40, e_rate40, moment_diff40, b_value40, d_value40];
i45=[cp45, mc_rate45, mc_var45, shear_fraction45, e_var45, e_rate45, moment_diff45, b_value45, d_value45];
i50=[cp50, mc_rate50, mc_var50, shear_fraction50, e_var50, e_rate50, moment_diff50, b_value50, d_value50];

input_master1=[i0; i5; i10; i15; i20; i30; i35; i40; i45; i50];
% input_master2=[i0(:,2:end); i5(:,2:end); i10(:,2:end); i15(:,2:end); i20(:,2:end); i30(:,2:end); i35(:,2:end); i40(:,2:end); i45(:,2:end); i50(:,2:end)];
% input_master3=[i0(:,[1,2,4,6,7]); i5(:,[1,2,4,6,7]); i10(:,[1,2,4,6,7]); i15(:,[1,2,4,6,7]); i20(:,[1,2,4,6,7]); i30(:,[1,2,4,6,7]); i35(:,[1,2,4,6,7]); i40(:,[1,2,4,6,7]); i45(:,[1,2,4,6,7]); i50(:,[1,2,4,6,7])];
% input_master4=[i0(:,[2,4,6,7]); i5(:,[2,4,6,7]); i10(:,[2,4,6,7]); i15(:,[2,4,6,7]); i20(:,[2,4,6,7]); i30(:,[2,4,6,7]); i35(:,[2,4,6,7]); i40(:,[2,4,6,7]); i45(:,[2,4,6,7]); i50(:,[2,4,6,7])];
%%
%output arrays
% o0=[stress_tf5, strain_tf5, cp5];
% o5=[stress_tf5, strain_tf5, cp5];
% o10=[stress_tf10, strain_tf10, cp10];
% o15=[stress_tf15, strain_tf15, cp15];
% o20=[stress_tf20, strain_tf20, cp20];
% o30=[stress_tf30, strain_tf30, cp30];
% o35=[stress_tf35, strain_tf35, cp35];
% o40=[stress_tf40, strain_tf40, cp40];
% o45=[stress_tf45, strain_tf45, cp45];
% o50=[stress_tf50, strain_tf50, cp50];
o_stress_tf=[stress_tf0; stress_tf5; stress_tf10; stress_tf15; stress_tf20; stress_tf30; stress_tf35; stress_tf40; stress_tf45; stress_tf50];
o_strain_tf=[strain_tf0; strain_tf5; strain_tf10; strain_tf15; strain_tf20; strain_tf30; strain_tf35; strain_tf40; strain_tf45; strain_tf50];
o_cp=[cp0; cp5; cp10; cp15; cp20; cp30; cp35; cp40; cp45; cp50];
%%
%25 mpa - testing set inputs
i25=[cp25, mc_rate25, mc_var25, shear_fraction25, e_var25, e_kur25, moment_diff25, b_value25, d_value25];
input_25t_master1=i25;
% input_25t_master2=i25(:,(2:end));
% input_25t_master3=i25(:,[1,2,4,6,7]);
% input_25t_master4=i25(:,[2,4,6,7]);
%25 mpa - testing set outputs
o_25t_stress_tf=stress_tf25;
o_25t_strain_tf=strain_tf25;
o_25t_cp=cp25;
%%
run p3_prediction_RF


