#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Container/Vector.h"
#include "/home/codeleaded/System/Static/Library/Victor2.h"
#include "/home/codeleaded/System/Static/Library/Vdctor2.h"
#include "/home/codeleaded/System/Static/Library/Complex2.h"
#include "/home/codeleaded/System/Static/Library/TransformedView.h"
#include "/home/codeleaded/System/Static/Library/Thread.h"
#include "/home/codeleaded/System/Static/Library/ThreadPool.h"

#define nMaxThreads 32

TransformedViewD tv;
int* pFractal = NULL;
int nMode = 0;
int nIterations = 1;

void CreateFractalBasic(const Vic2* pix_tl, const Vic2* pix_br, const Vdc2* frac_tl, const Vdc2* frac_br, const int iterations){
	double x_scale = (frac_br->x - frac_tl->x) / ((double)pix_br->x - (double)pix_tl->x);
	double y_scale = (frac_br->y - frac_tl->y) / ((double)pix_br->y - (double)pix_tl->y);
	
	for (int y = pix_tl->y; y < pix_br->y; y++)
	{
		for (int x = pix_tl->x; x < pix_br->x; x++)
		{
			Complex2 c = { x * x_scale + frac_tl->x, y * y_scale + frac_tl->y };
			Complex2 z = { x * x_scale + frac_tl->x, y * y_scale + frac_tl->y };
			int n = 0;
			while (Complex2_Mag(z) < 2.0 && n < iterations)
			{
				z = Complex2_Add(Complex2_Mul(z,z),c);
				n++;
			}
			pFractal[y * GetWidth() + x] = n;
		}
	}
}
void CreateFractalPreCalculate(const Vic2* pix_tl, const Vic2* pix_br, const Vdc2* frac_tl, const Vdc2* frac_br, const int iterations){
	double x_scale = (frac_br->x - frac_tl->x) / ((double)pix_br->x - (double)pix_tl->x);
	double y_scale = (frac_br->y - frac_tl->y) / ((double)pix_br->y - (double)pix_tl->y);

	double x_pos = frac_tl->x;
	double y_pos = frac_tl->y;

	int y_offset = 0;
	int row_size = pix_br->x - pix_tl->x;

	int x, y, n;
	Complex2 c = { 0.0f,0.0f };
	Complex2 z = { 0.0f,0.0f };

	for (y = pix_tl->y; y < pix_br->y; y++)
	{
		x_pos = frac_tl->x;
		for (x = pix_tl->x; x < pix_br->x; x++)
		{				
			c = (Complex2){ x_pos, y_pos };
			z = (Complex2){ 0,0 };

			n = 0;
			while (Complex2_Mag(z) < 2.0 && n < iterations)
			{
				z = Complex2_Add(Complex2_Mul(z,z),c);
				n++;
			}

			pFractal[y_offset + x] = n;
			x_pos += x_scale;
		}

		y_pos += y_scale;
		y_offset += row_size;
	}
}
void CreateFractalNoComplex2(const Vic2* pix_tl, const Vic2* pix_br, const Vdc2* frac_tl, const Vdc2* frac_br, const int iterations){
	double x_scale = (frac_br->x - frac_tl->x) / ((double)pix_br->x - (double)pix_tl->x);
	double y_scale = (frac_br->y - frac_tl->y) / ((double)pix_br->y - (double)pix_tl->y);
	double x_pos = frac_tl->x;
	double y_pos = frac_tl->y;
	int y_offset = 0;
	int row_size = GetWidth();
	int x, y, n;
	double cr = 0;
	double ci = 0;
	double zr = 0;
	double zi = 0;
	double re = 0;
	double im = 0;
	for (y = pix_tl->y; y < pix_br->y; y++)
	{
		x_pos = frac_tl->x;
		ci = y_pos;
		for (x = pix_tl->x; x < pix_br->x; x++)
		{
			cr = x_pos;
			zr = 0;
			zi = 0;
			n = 0;
			while ((zr * zr + zi * zi) < 4.0 && n < iterations)
			{
				re = zr * zr - zi * zi + cr;
				im = zr * zi * 2.0 + ci;
				zr = re;
				zi = im;
				n++;
			}
			pFractal[y_offset + x] = n;
			x_pos += x_scale;
		}
		y_pos += y_scale;
		y_offset += row_size;
	}
}
void CreateFractalIntrinsics(const Vic2* pix_tl, const Vic2* pix_br, const Vdc2* frac_tl, const Vdc2* frac_br, const int iterations){
	double x_scale = (frac_br->x - frac_tl->x) / ((double)pix_br->x - (double)pix_tl->x);
	double y_scale = (frac_br->y - frac_tl->y) / ((double)pix_br->y - (double)pix_tl->y);
	double y_pos = frac_tl->y;
	int y_offset = 0;
	int row_size = GetWidth();
	int x, y;
	__m256d _a, _b, _two, _four, _mask1;
	__m256d _zr, _zi, _zr2, _zi2, _cr, _ci;
	__m256d _x_pos_offsets, _x_pos, _x_scale, _x_jump;
	__m256i _one, _c, _n, _iterations, _mask2;
	_one = _mm256_set1_epi64x(1);
	_two = _mm256_set1_pd(2.0);
	_four = _mm256_set1_pd(4.0);
	_iterations = _mm256_set1_epi64x(iterations);
	_x_scale = _mm256_set1_pd(x_scale);
	_x_jump = _mm256_set1_pd(x_scale * 4);
	_x_pos_offsets = _mm256_set_pd(0, 1, 2, 3);
	_x_pos_offsets = _mm256_mul_pd(_x_pos_offsets, _x_scale);
	
	for (y = pix_tl->y; y < pix_br->y; y++){
		_a = _mm256_set1_pd(frac_tl->x);
		_x_pos = _mm256_add_pd(_a, _x_pos_offsets);
		_ci = _mm256_set1_pd(y_pos);
		
		for (x = pix_tl->x; x < pix_br->x; x += 4){
			_cr = _x_pos;
			_zr = _mm256_setzero_pd();
			_zi = _mm256_setzero_pd();
			_n = _mm256_setzero_si256();
			
			do{
				_zr2 = _mm256_mul_pd(_zr, _zr);
				_zi2 = _mm256_mul_pd(_zi, _zi);
				_a = _mm256_sub_pd(_zr2, _zi2);
				_a = _mm256_add_pd(_a, _cr);
				_b = _mm256_mul_pd(_zr, _zi);
				//_b = _mm256_fmadd_pd(_b, _two, _ci);
				_b = _mm256_mul_pd(_b, _two);
				_b = _mm256_add_pd(_b,_ci);
				_zr = _a;
				_zi = _b;
				_a = _mm256_add_pd(_zr2, _zi2);
				_mask1 = _mm256_cmp_pd(_a, _four, _CMP_LT_OQ);
				_mask2 = _mm256_cmpgt_epi64(_iterations, _n);
				_mask2 = _mm256_and_si256(_mask2, _mm256_castpd_si256(_mask1));
				_c = _mm256_and_si256(_one, _mask2);											
				_n = _mm256_add_epi64(_n,_c);
			}while(_mm256_movemask_pd(_mm256_castsi256_pd(_mask2)) > 0);

			pFractal[y_offset + x + 0] = (int)_mm256_extract_epi64(_n,3);
			pFractal[y_offset + x + 1] = (int)_mm256_extract_epi64(_n,2);
			pFractal[y_offset + x + 2] = (int)_mm256_extract_epi64(_n,1);
			pFractal[y_offset + x + 3] = (int)_mm256_extract_epi64(_n,0);
			_x_pos = _mm256_add_pd(_x_pos, _x_jump);
		}
		y_pos += y_scale;
		y_offset += row_size;
	}
}



typedef struct Element{
	Thread hThread;
	char bLocked;
	char bAlive;
	Vic2 a;
	Vic2 b;
	Vdc2 c;
	Vdc2 d;
	int i;
} Element;

void* Thread_Func(void* arg){
	Element* e = arg;
	CreateFractalIntrinsics(&e->a,&e->b,&e->c,&e->d,e->i);
	return NULL;
}
void CreateFractalThreads(const Vic2* pix_tl, const Vic2* pix_br, const Vdc2* frac_tl, const Vdc2* frac_br, const int iterations){
	int nSectionWidth = (pix_br->x - pix_tl->x) / nMaxThreads;
	double dFractalWidth = (frac_br->x - frac_tl->x) / (double)nMaxThreads;
	
	Element Args[nMaxThreads];
	for (size_t i = 0; i < nMaxThreads; i++){
		Args[i].a = (Vic2){ pix_tl->x + nSectionWidth * (i), pix_tl->y };
		Args[i].b = (Vic2){ pix_tl->x + nSectionWidth * (i + 1), pix_br->y };
		Args[i].c = (Vdc2){ frac_tl->x + dFractalWidth * (double)(i), frac_tl->y };
		Args[i].d = (Vdc2){ frac_tl->x + dFractalWidth * (double)(i + 1), frac_br->y };
		Args[i].i = (int)iterations;
		Args[i].hThread = Thread_New(NULL,Thread_Func,&Args[i]);
		Thread_Start(&Args[i].hThread);
	}
	
	for (size_t i = 0; i < nMaxThreads; i++){
		Thread_Join(&Args[i].hThread);
	}
}


typedef struct PoolElement{
	Vic2 a;
	Vic2 b;
	Vdc2 c;
	Vdc2 d;
	int i;
} PoolElement;

ThreadPool pool;

void* Thread_Pool_Func(ThreadWorker* tw){
	PoolElement* e = (PoolElement*)tw->data;
	CreateFractalIntrinsics(&e->a,&e->b,&e->c,&e->d,e->i);
	return NULL;
}
void CreateFractalThreadPool(const Vic2* pix_tl, const Vic2* pix_br, const Vdc2* frac_tl, const Vdc2* frac_br, const int iterations){
	const int nSectionWidth = (pix_br->x - pix_tl->x) / pool.size;
	const double dFractalWidth = (frac_br->x - frac_tl->x) / (double)pool.size;
	for(int i = 0;i < pool.size;i++){
		ThreadWorker* tw = (ThreadWorker*)Vector_Get(&pool,i);
		PoolElement* pe = (PoolElement*)tw->data;
		pe->a = (Vic2){ pix_tl->x + nSectionWidth * (i), pix_tl->y };
		pe->b = (Vic2){ pix_tl->x + nSectionWidth * (i + 1), pix_br->y };
		pe->c = (Vdc2){ frac_tl->x + dFractalWidth * (double)(i), frac_tl->y };
		pe->d = (Vdc2){ frac_tl->x + dFractalWidth * (double)(i + 1), frac_br->y };
		pe->i = (int)iterations;
	}

	// const int nSectionHeight = (pix_br->y - pix_tl->y) / pool.size;
	// const double dFractalHeight = (frac_br->y - frac_tl->y) / (double)pool.size;
	// for(int i = 0;i < pool.size;i++){
	// 	ThreadWorker* tw = (ThreadWorker*)Vector_Get(&pool,i);
	// 	PoolElement* pe = (PoolElement*)tw->data;
	// 	pe->a = (Vic2){ pix_tl->x,pix_tl->y + nSectionHeight * (i) };
	// 	pe->b = (Vic2){ pix_br->x,pix_tl->y + nSectionHeight * (i + 1) };
	// 	pe->c = (Vdc2){ frac_tl->x,frac_tl->y + dFractalHeight * (double)(i) };
	// 	pe->d = (Vdc2){ frac_br->x,frac_tl->y + dFractalHeight * (double)(i + 1) };
	// 	pe->i = (int)iterations;
	// }

	ThreadPool_Start(&pool);
	ThreadPool_Join(&pool);
}


void Setup(AlxWindow* w){
	tv = TransformedViewD_New((Vdc2){ GetWidth(),GetHeight() });
	
	pFractal = (int*)malloc((size_t)GetWidth() * (size_t)GetHeight() * sizeof(int));
	memset(pFractal,0,(size_t)GetWidth() * (size_t)GetHeight() * sizeof(int));

	pool = ThreadPool_New((void*(*)(void*))Thread_Pool_Func,sizeof(PoolElement),nMaxThreads);

	ResizeAlxFont(16,16);
}
void Update(AlxWindow* w){
	tv.ZoomSpeed = (float)w->ElapsedTime * 5.0f;
	TransformedViewD_HandlePanZoom(&tv,window.Strokes,(Vdc2){ GetMouse().x,GetMouse().y });
	
	Vic2 pix_tl = { 0,0 };
	Vic2 pix_br = { GetWidth(), GetHeight() };
	Vdc2 frac_tl = { -2, -1 };
	Vdc2 frac_br = { 1, 1 };
	
	frac_tl = TransformedViewD_ScreenWorldPos(&tv,(Vdc2){ pix_tl.x,pix_tl.y });
	frac_br = TransformedViewD_ScreenWorldPos(&tv,(Vdc2){ pix_br.x,pix_br.y });

	if (Stroke(ALX_KEY_1).PRESSED) 	nMode = 0;
	if (Stroke(ALX_KEY_2).PRESSED) 	nMode = 1;
	if (Stroke(ALX_KEY_3).PRESSED) 	nMode = 2;
	if (Stroke(ALX_KEY_4).PRESSED) 	nMode = 3;
	if (Stroke(ALX_KEY_5).PRESSED) 	nMode = 4;
	if (Stroke(ALX_KEY_6).PRESSED) 	nMode = 5;
	if (Stroke(ALX_KEY_UP).DOWN) 	nIterations += 1;
	if (Stroke(ALX_KEY_DOWN).DOWN) 	nIterations -= 1;
	if (nIterations < 1) nIterations = 1;
	
	Timepoint tp1 = Time_Nano();
	
	switch (nMode){
		case 0: CreateFractalBasic(&pix_tl,&pix_br,&frac_tl,&frac_br,nIterations); 			break;
		case 1: CreateFractalPreCalculate(&pix_tl,&pix_br,&frac_tl,&frac_br,nIterations); 	break;
		case 2: CreateFractalNoComplex2(&pix_tl,&pix_br,&frac_tl,&frac_br,nIterations); 	break;
		case 3: CreateFractalIntrinsics(&pix_tl,&pix_br,&frac_tl,&frac_br,nIterations); 	break;
		case 4: CreateFractalThreads(&pix_tl,&pix_br,&frac_tl,&frac_br,nIterations); 		break;
		case 5: CreateFractalThreadPool(&pix_tl,&pix_br,&frac_tl,&frac_br,nIterations); 	break;
	}
	
	Timepoint tp2 = Time_Nano();
	double elapsedTime = Time_Elapsed(tp1,tp2);

	for (int i = 0; i < GetWidth() * GetHeight(); i++){
		int it = pFractal[i];
		float n = (float)it;
		float a = 0.05f;

		window.Buffer[i] = Pixel_toRGBA(
			0.5f * F32_Sin_Cb(a * n) + 0.5f,
			0.5f * F32_Sin_Cb(a * n + 2.094f) + 0.5f,
			0.5f * F32_Sin_Cb(a * n + 4.188f) + 0.5f,
			1.0f
		);
	}

	switch (nMode){
		case 0: RenderCStr("1) Naive Method",					0,0,WHITE); break;
		case 1: RenderCStr("2) Precalculate Method",			0,0,WHITE); break;
		case 2: RenderCStr("3) Hand-code Maths Method",			0,0,WHITE); break;
		case 3: RenderCStr("4) Vector Extensions (AVX2) Method",0,0,WHITE); break;
		case 4: RenderCStr("5) Threads Method",					0,0,WHITE); break;
		case 5: RenderCStr("6) ThreadPool Method",				0,0,WHITE); break;
	}

	String str = String_Format("| ES Taken: %fs - Iterations: %d |",elapsedTime,nIterations);
	char* cstr = String_CStr(&str);
	RenderCStr(cstr,0,30,WHITE);
	free(cstr);
	String_Free(&str);
}
void Delete(AlxWindow* w){
	ThreadPool_Free(&pool);
	free(pFractal);
}
void Resize(AlxWindow* w){
    if(pFractal) free(pFractal);
	pFractal = (int*)malloc((size_t)GetWidth() * (size_t)GetHeight() * sizeof(int));
}

int main(){
    if(CreateX("Game Test",1280,720,1,1,Setup,Update,Delete,Resize))
        Start();
    return 0;
}