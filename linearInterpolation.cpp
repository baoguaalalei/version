void linearInterpolation(int n_data,int *data,float *pWave_lineInter,int * pos_peak,int n_peak,int nSamplesPerSecond)
{
	if (n_peak<10)
	{
		free(pos_peak);
		return;
	}	
    double a,b;//一阶线性参数
    long x1,x2;//保存线性函数的两个已知x值
	long x;
	int pos_N;
	pos_N=1;//从第1个峰值位置开始
    while (pos_N<n_peak)
    {
		x1=pos_peak[pos_N-1];
		x2=pos_peak[pos_N];
		a=(double)(data[x1]-data[x2])/(x1-x2);
	    b=(x1-(double)data[x1]/data[x2]*x2)/(x1-x2);//数据比较大，超出了范围，此处采取的方法是先除以data[x2]再乘以该数的方法
		b*=data[x2];
		//printf("%f %f\n",a,b);
		x=x1;
		while (x<=x2)
		{
			//pWave_lineInter[x]=x*a+b;
			pWave_lineInter[x]=data[x]-(x*a+b);
			//printf("%f \n",pWave_lineInter[x]);
			x++;
		}
		pos_N++;
	}
	#if 1	
	FILE *fp;
    fp=fopen("D:\\实验结果3\\linearInterpolation_result.csv","w");
	int num=0;
	if(fp)
	{
		while (num<pos_peak[n_peak-1])
		{
			if (num<pos_peak[0])//没有参与线性差值的点直接赋予零值
			{
				fprintf(fp,"%f\n",0);
			}
			else
		    {
				fprintf(fp,"%f\n",pWave_lineInter[num]);
			}
			num++;
		}
		fclose(fp);		
	}
#endif
return;
};
