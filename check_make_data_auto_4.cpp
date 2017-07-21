void  check_make_data_auto(int n_data,int*data,int n_peak,int *pos_peak,int n_valley,
		int *pos_valley,int n_grade,int *pos_grade,int &nt,float * f_NN,float * f_NN_t)
{
	int i;
	double ratio_ave=0.0,ratio_var=0.0,ave_ave=0.0,ave_var=0.0;
	double ratio_ave_peak=0.0,ratio_var_peak=0.0,ave_ave_peak=0.0,ave_var_peak=0.0;
	double ratio_ave_grade=0.0,ratio_var_grade=0.0,ave_ave_grade=0.0,ave_var_grade=0.0;
	double ratio_ave_valley=0.0,ratio_var_valley=0.0,ave_ave_valley=0.0,ave_var_valley=0.0;

	// 数据点太少，无法处理
	if (n_valley<10) return;
	// 计算谷值的均值，方差	
	/************************************************************************/
	/* 用于检测间期的变化是否符合规律正态分布                                                                     */
	/************************************************************************/
 
	for (i=1;i<n_peak-1;i++)
	{
		double ave_peak = pos_peak[i+1]-pos_peak[i-1];
		double ratio_peak= (double)(pos_peak[i+1]-pos_peak[i])/(pos_peak[i]-pos_peak[i-1]);
		ave_ave_peak+=ave_peak;
		ave_var_peak+=ave_peak*ave_peak;
		ratio_ave+=ratio_peak;
		ratio_var_peak+=ratio_peak*ratio_peak;
	}
	ave_ave_peak/=(n_peak-2);
	ave_var_peak=sqrt(ave_var_peak/(n_peak-2)-ave_ave_peak*ave_ave_peak);
	ratio_ave_peak/=(n_peak-2);		// 数据少的时候应该用中值，均值很容易受干扰!!! 当间隔为1的错误时
	ratio_var_peak=sqrt(ratio_var_peak/(n_peak-2)-ratio_ave_peak*ratio_ave_peak);
	printf("峰值间期均值=%f,峰值间期均方差=%f,峰值比值均值=%f,峰值比值均方差=%f\n",(float)ave_ave_peak,(float)ave_var_peak,(float)ratio_ave_peak,(float)ratio_var_peak);
	
	for (i=1;i<n_grade-1;i++)
	{
		double ave_grade = pos_grade[i+1]-pos_grade[i-1];
		double ratio_grade = (double)(pos_grade[i+1]-pos_grade[i])/(pos_grade[i]-pos_grade[i-1]);
		ave_ave_grade+=ave_grade;
		ave_var_grade+=ave_grade*ave_grade;
		ratio_ave_grade+=ratio_grade;
		ratio_var_grade+=ratio_grade*ratio_grade;
	}
	ave_ave_grade/=(n_grade-2);
	ave_var_grade=sqrt(ave_var_grade/(n_grade-2)-ave_ave_grade*ave_ave_grade);
	ratio_ave_grade/=(n_grade-2);		// 数据少的时候应该用中值，均值很容易受干扰!!! 当间隔为1的错误时
	ratio_var_grade=sqrt(ratio_var_grade/(n_grade-2)-ratio_ave_grade*ratio_ave_grade);
	printf("梯度间期均值=%f,梯度间期均方差=%f,梯度比值均值=%f,梯度比值均方差=%f\n",(float)ave_ave_grade,(float)ave_var_grade,(float)ratio_ave_grade,(float)ratio_var_grade);
	
	for(i=1;i<n_valley-1;i++)
	{
		double ave_valley=pos_valley[i+1]-pos_valley[i-1];
		double ratio_valley=(double)(pos_valley[i+1]-pos_valley[i])/(pos_valley[i]-pos_valley[i-1]);
		//		printf("%d,%d,%f\n",i,pos_grade[i+1]-pos_grade[i],float(ratio));
		ave_ave_valley+=ave_valley;      //ava_ava是间期均值，ave_var是间期均方差，
		ave_var_valley+=ave_valley*ave_valley;
		ratio_ave_valley+=ratio_valley;  //ratio_ave是比值均值，ratio_var是比值均方差
		ratio_var_valley+=ratio_valley*ratio_valley;
	}
	ave_ave_valley/=(n_valley-2);	//间期的平均值的两倍
	ave_var_valley=sqrt(ave_var_valley/(n_valley-2)-ave_ave_valley*ave_ave_valley);
	ratio_ave_valley/=(n_valley-2);		// 数据少的时候应该用中值，均值很容易受干扰!!! 当间隔为1的错误时
	ratio_var_valley=sqrt(ratio_var_valley/(n_valley-2)-ratio_ave_valley*ratio_ave_valley);
	printf("波谷间期均值=%f,波谷间期均方差=%f,波谷比值均值=%f,波谷比值均方差=%f\n",(float)ave_ave_valley,(float)ave_var_valley,(float)ratio_ave_valley,(float)ratio_var_valley);
	
	//比较波峰、波谷、最大梯度，均方差最小的相对稳定，错误的点比较少，作为筛选的条件

	ave_ave=ave_ave_grade;
	ave_var=ave_var_grade;
	if (ave_var>ave_var_valley)
	{
		ave_var=ave_var_valley;
		ave_ave=ave_ave_valley;
	}
	if (ave_var>ave_var_peak)
	{
		ave_var=ave_var_peak;
		ave_ave=ave_ave_peak;
	}
	

	//进行间期变化的统计分析，通过excel表格进行分析发现，数据分布符合正态分布
	double ratio_ave_RR=0.0,ratio_var_RR=0.0,ave_ave_RR=0.0,ave_var_RR=0.0;
	for (i=1;i<n_peak-1;i++)
	{
		double ave_RR=pos_peak[i+1]+pos_peak[i-1]-2*pos_peak[i];
	    ave_ave_RR+=ave_RR;
		ave_var_RR+=ave_RR*ave_RR;
	}
	ave_ave_RR/=(n_peak-2);
	ave_var_RR=sqrt(ave_var_RR/(n_peak-2)-ave_ave_RR*ave_ave_RR);	

	int p_grade, p_valley,p_peak;//用于保存当前位置
	int i_g;//最大梯度点寻找范围	
	int i_v;//波谷点寻找的范围
	int i_p;//峰值点寻找范围
	bool grade_find;//重新定位后第一个最大梯度点找到的标志
	bool peak_find;//重新定位后第一个峰值点找到的标志
	bool valley_find;//重新定位后第一个波谷点找到的标志
    long  t;//贯穿整个搜索过程的时间
	
	//寻找有效间期序列
	//从波谷点的n_valley=1开始找
	nt=0;
	t=pos_valley[1];
	
	//用于标注找到的波谷，最大梯度，波谷值是否是坏点
	bool vBad;
	bool gBad;
	bool pBad;	

	while(1)
    {
		//先找波谷，再找最大梯度，最后找波峰
		valley_find=false;
		for(i=0;i<n_valley-1;i++)
		{
			if(pos_valley[i]>t&&(pos_valley[i]-t)>ave_ave/4)
			{
				valley_find=true;
				p_valley=i;
				break;		
			}
		}
		if (!valley_find)//没有找到一个符合要求的波谷值，直接跳出循环
		{
			break;
		}
		vBad=true;
		vBad=check_bad_data(p_valley,pos_valley,ave_ave,ave_var);//检查点是否有问题
		if (!vBad)//符合要求
		{
			f_NN[nt]=pos_valley[p_valley]-pos_valley[p_valley-1];
			t=pos_valley[p_valley];
			printf("%d %f %d valley\n",nt,f_NN[nt],t);
			nt++;
		} 
		else
		{
			grade_find=false;
			for(i=0;i+1<n_grade-1;i++)
			{
				if(pos_grade[i]>t)
				{
					grade_find=true;
					p_grade=i;
					break;		
				}
			}
			if (!grade_find)//没有找到大于标度位置的最大梯度值，直接跳出循环，如果可以找到，那么找到后再做进一步的判断
			{
				break;
			}
			gBad=true;
			grade_find=(pos_grade[p_grade+1]-t)<ave_ave;
			if (grade_find)
			{
				gBad=check_bad_data(p_grade+1,pos_grade,ave_ave,ave_var);//检查点是否有问题.这里需要注意最大梯度点的范围，p_grade是否越界
			}
		
		    if (!gBad)
			{
				f_NN[nt]=pos_grade[p_grade+1]-pos_grade[p_grade];
				t+=f_NN[nt];
				printf("%d %f %d grade\n",nt,f_NN[nt],t);
				nt++;
			} 
			else
			{//最大梯度不符合要求时采用峰峰值来判断
				peak_find=false;
				for (i=1;i+1<n_peak-1;i++)
				{
					if (pos_peak[i]>t)
					{
						peak_find=true;
						p_peak=i;
						break;
					}
				}
				if (!peak_find)
				{
					break;
				}
				pBad=true;
				peak_find=(pos_peak[p_peak+1]-t)>ave_ave/2;//后期会用到该值
				if(peak_find)
				{
					pBad=check_bad_data(p_peak+1,pos_peak,ave_ave,ave_var);				
				}
				
				if (!pBad)
				{
					f_NN[nt]=pos_peak[p_peak+1]-pos_peak[p_peak];
					t+=f_NN[nt];
					printf("%d %f %d peak\n",nt,f_NN[nt],t);
					nt++;
				}
				else//当波谷、最大梯度、波峰值均不符合要求时，寻找第一个大于t的valley重新给t赋值，如果没有找到合适的，跳出循环
				{
					//最后需要判断一下是不是被勿删掉的数据，根据间期的变化是否符合规律来判断，按照波谷、最大梯度、波峰依次做判断，只要有两个可以解决这个
					int count_RR=0;//标记波谷、波峰、最大梯度三个值中符合要求的点的个数
					vBad=true;
					gBad=true;
					pBad=true;	
					int RR1=0,RR2=0,RR3=0;
					int RR=0;
					if (p_grade<n_grade-3&&p_grade>=1)//间期变化规律，n_grade-3是最大梯度位置最后一个d(RR)
					{
						if (grade_find)
						{
							gBad=check_bad_RR(p_grade+1,pos_grade,ave_ave_RR,ave_var_RR);
						}
						if (!gBad)
						{
							count_RR++;
							RR1=pos_grade[p_grade+1]-pos_grade[p_grade];
						}						
					}
					if (p_valley<n_valley-2&&p_valley>=2)
					{
						if (valley_find)
						{
							vBad=check_bad_RR(p_valley,pos_valley,ave_ave_RR, ave_var_RR);
						}
						if (!vBad)
						{
							count_RR++;
							RR2=pos_valley[p_valley]-pos_valley[p_valley-1];
						}
					}
					if (p_peak<n_peak-3&&p_peak>=1)
					{
						if (peak_find)
						{
							pBad=check_bad_RR(p_peak+1,pos_peak,ave_ave_RR,ave_var_RR);							
						}
						if (!pBad)
						{
							count_RR++;
							RR3=pos_peak[p_peak+1]-pos_peak[p_peak];
						}
					}
					if (count_RR>=2)
					{					
						if (abs(RR1-RR2)<6)//间期变化的差异选择哪个值，需要进一步推敲，暂定为6,这样的心跳误差大约为0.42次/min
						{
							RR=(RR1+RR2)/2;
						}
						else if (abs(RR1-RR3)<6)
						{
							RR=(RR1+RR3)/2;
						}
						else if (abs(RR2-RR3)<6)
						{
							RR=(RR2+RR3)/2;
						}
					}
					if(RR>0)
					{
						f_NN[nt]=RR;
						t+=f_NN[nt];
						printf("%d %f %d peak\n",nt,f_NN[nt],t);
						nt++;
					}
					else
					{
						valley_find=false;
						for (i=1;i<n_valley-1;i++)
						{
							if (pos_valley[i]>t)
							{
								valley_find=true;
								p_valley=i;
								break;
							}												
						}
						if (!valley_find)
						{
							break;
						}
						else
						{
							t=pos_valley[p_valley];
						}
					
					}				
				}
			}			
		}
	}	
	f_NN_t[0]=pos_valley[1];
	for (i=0;i<nt;i++)
	{
		f_NN_t[i+1]=f_NN[i]+f_NN_t[i];
	}
	
	FILE *fp_t;
	fp_t=fopen("D:\\me_2017-04-17 14_22_55_valley_result_ok.csv","w");
	long num=0;
	if(fp_t)
	{
		while (num<=nt)
		{
			fprintf(fp_t,"%f\n",f_NN_t[num]);
			num++;
		}
	}
}
