int PulseHarmonic(float *pulse,int n_pulse,float fSampleRate, float *pulseOut,int n_FFT){
	float sum;
	int i;	
	float *t =(float*)malloc(sizeof(float)*n_pulse);            //原pulse的x轴		   
	float *pulse_resample_y=(float*)malloc(sizeof(float)*n_FFT);//重采样的y值
	float *pulse_resample_x=(float*)malloc(sizeof(float)*n_FFT);//重采样的x值
	float fs_new=n_FFT*fSampleRate/n_pulse;                     //重采样的频率步长，也就是脉搏信号第一谐波频率
	//检验内存分配
	if (NULL==t||NULL==pulse_resample_y||NULL==pulse_resample_x)
	{
		printf("Memory failed!");
		if (t) free(t);
		if (pulse_resample_x) free(pulse_resample_x);
		if (pulse_resample_y) free(pulse_resample_y);
		return -1;
	}	
	sum=0.0;
	for (i=0;i<n_pulse;i++)	
		sum+=pulse[i];
	sum/=n_pulse;
	for (i=0;i<n_pulse;i++) 
		pulse[i]-=sum;	
	for (i=0;i<n_pulse;i++) t[i]=i;
	float dx=(float)(n_pulse-1)/(n_FFT-1);                      //重采样单步长
	for (int i=0;i<n_FFT;i++)                                   //输出时间x的初始化
		pulse_resample_x[i]=dx*i;	

	if (Spline_xwh(t,pulse,n_pulse,pulse_resample_x,pulse_resample_y,n_FFT))
	{
		free(t);
		free(pulse_resample_x);
		free(pulse_resample_y);
		return -1;
	}
	//FFT变换	
	fft_main(pulse_resample_y,n_FFT,pulseOut,n_FFT,fs_new,0,NULL,1); 
	for (i=0;i<n_FFT;i++)
		pulseOut[i]=sqrt(pulseOut[i])/norm;	

	free(t);
    free(pulse_resample_x);
	free(pulse_resample_y);
	return 0;
}
