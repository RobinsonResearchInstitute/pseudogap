//Writen by James Storey, Robinson Research Institute, 2020.

let previousProgress = 0;

const kb_ = 8.62e-5;
const alphaBCS_ = 2.14;

let alpha_ = alphaBCS_;
let Tc_ = 80.0;
let Eg_ = 2.0*alphaBCS_*kb_*Tc_;
let thc_ = Math.PI/4.0;



let energy_ = [];
let dos_ = [];
let NSdos_=[];
let E1_ = -0.1;
let E2_ = 0.1;
const NE_ = 201;
let dE_ = (E2_-E1_)/(NE_-1);

let gap_ = [];
let entropy_ = [];
let NSentropy_ = [];
let susceptibility_ = [];
let NSsusceptibility_ = [];
let temperature_ = [];
let gamma_ = [];
let NSgamma_ = [];
let DF_ = [];
let Dt_ = [];
let lambda_ = [];
let T1_ = 0;
//let T2_ = 30;
const NT_ = 301;
let nTc_ = 80;//0.8*(NT_-1);
let dT_ = (Tc_-T1_)/nTc_;
let T2_ = T1_+(NT_-1)*dT_;




onmessage = function(event) {
    // Perform the prime number search.
    Tc_ = event.data.Tc;
    this.console.log(Tc_);
    alpha_ = event.data.alpha;
    Eg_ = event.data.Eg*alphaBCS_*kb_*Tc_;
    thc_ = event.data.thc*Math.PI/180.0;
    initArrays();
    computeDos();
    postMessage(
        {messageType: "Dos", E: energy_, N: dos_, Nn: NSdos_}
      );

    calculateDelta();  
    postMessage(
        {messageType: "Delta", D: gap_}
      );
    calculateEntropy(alpha_, Tc_);  // Send back the results.
    postMessage(
        {messageType: "Entropy", S: entropy_, Sn: NSentropy_}
      );
    calculateGamma();
    postMessage(
        {messageType: "Gamma", g: gamma_, gn: NSgamma_}
      );
    calculateSusceptibility(alpha_, Tc_);
    postMessage(
        {messageType: "Susceptibility", X: susceptibility_, Xn: NSsusceptibility_}
      );
    calculateDF();
    postMessage(
        {messageType: "DF", F: DF_}
      );
    calculateDt();
    postMessage(
        {messageType: "Dt", Dt: Dt_}
      );
    calculateLvsT();
    postMessage(
        {messageType: "Lambda", L: lambda_}
      );  
    
  };
  

  function initArrays(){
    energy_ = new Array(NE_);
    for(let i = 0; i<NE_; i++){
        energy_[i] = (E1_*100+i*dE_*100)/100;
    }    
    dos_ = new Array(NE_).fill(0);
    NSdos_ = new Array(NE_).fill(0);

    temperature_ = new Array(NT_);
    for(let i = 0; i<NT_; i++){
        temperature_[i] = T1_+i*dT_;
    }
    gap_ = new Array(NT_).fill(0);
    entropy_ = new Array(NT_).fill(0);
    NSentropy_ = new Array(NT_).fill(0);
    susceptibility_ = new Array(NT_).fill(0);
    NSsusceptibility_ = new Array(NT_).fill(0);
    //gamma_ = new Array(NT_).fill(0);
    //NSgamma_ = new Array(NT_).fill(0);
    DF_ = new Array(NT_).fill(0);
    Dt_ = new Array(NT_).fill(0);
    lambda_ = new Array(NT_).fill(0);
}


function computeDos(){
	
    let th1 = 0;
    let th2 = Math.PI/4.0;
    const Nth = 200;
    let dth = (th2-th1)/(Nth-1);
    let th = 0;
    var E = E1_;
    var N = 2.0;
    var delta = alpha_*kb_*Tc_;
    var deltatheta = 0;
    var N2 = 0;

    var x;
    for(let i = 0; i<NE_; i++){
        E = energy_[i];
        N = 0;
        N2 = 0;
        th = th1;
        deltatheta = delta*Math.cos(2.0*th);
        x = Math.pow(E,2)-Math.pow(deltatheta,2)-Math.pow(Eg_*Math.cos(2.0*th),2);
        if(x>0) N += 0.5*Math.abs(E)/Math.sqrt(x);
        if(E==0 && x==0){
                N+= 0.5;
        }
        x = Math.pow(E,2)-Math.pow(Eg_*Math.cos(2.0*th),2);
        if(x>0) N2 += 0.5*Math.abs(E)/Math.sqrt(x);
        if(E==0 && x==0){
                N2+=0.5;
            }

        th = th1+(Nth-1)*dth;
        deltatheta = delta*Math.cos(2.0*th);
        x = Math.pow(E,2)-Math.pow(deltatheta,2)-Math.pow(Eg_*Math.cos(2.0*th),2);
        if(x>0) N += 0.5*Math.abs(E)/Math.sqrt(x);
        if(E==0 && x==0){
                N+= 0.5;
        }
        x = Math.pow(E,2)-Math.pow(Eg_*Math.cos(2.0*th),2);
        if(x>0) N2 += 0.5*Math.abs(E)/Math.sqrt(x);
        if(E==0 && x==0){
                N2+=0.5;
            }

        for(let j = 1; j<Nth-1; j++){
            th = th1+j*dth;
            deltatheta = delta*Math.cos(2.0*th);
            x = Math.pow(E,2)-Math.pow(deltatheta,2)-Math.pow(Eg_*Math.cos(2.0*th),2);
            if(x>0) N += Math.abs(E)/Math.sqrt(x);
            if(E==0 && x==0){
                N+= 1.0;
            }
            x = Math.pow(E,2)-Math.pow(Eg_*Math.cos(2.0*th),2);
            if(x>0) N2 += Math.abs(E)/Math.sqrt(x);
            //N = fermi(E,0,0);
            if(E==0 && x==0){
                N2+=1.0;
            }
        }
        dos_[i]=(N*dth*4.0/Math.PI);
        NSdos_[i]=(N2*dth*4.0/Math.PI);
    }
    
}



  function calculateGamma(){
    for(let i = 0; i<NT_-1; i++){
        gamma_[i]=(entropy_[i+1]-entropy_[i])/dT_;
        NSgamma_[i]=(NSentropy_[i+1]-NSentropy_[i])/dT_;
    }
}

function calculateDt(){
    let x = 0;
    for(let i = 0; i<NT_; i++){
        x = 0;
        if(temperature_[i]<Tc_)
        x = Math.sqrt(DF_[i]/DF_[0])-(1.0-Math.pow(temperature_[i]/Tc_,2));
        Dt_[i] = x;
    }
}

function calculateDF(){
    let sum = 0;
    for(let i = nTc_-1; i>=0; i--){
        sum+=0.5*(entropy_[i]-NSentropy_[i]+entropy_[i+1]-NSentropy_[i+1]);
        DF_[i]=-sum*dT_;
    }
}



  function calculateEntropy(alpha, Tc) {
    // (The boring prime number calculations go in this function.)
    // Calculate the progress percentage.
    let progress = 0;//Math.round(i/list.length*100);  // Only send a progress update if the progress has changed
    // at least 1%.
    for(let i = 0; i<NT_; i++){
        progress = Math.floor(i/(NT_-1)*100);
        if(temperature_[i]==0){
            entropy_[i] = 0;
            NSentropy_[i] = 0;
        }else{
            entropy_[i] = S(alpha,Tc,temperature_[i]);
            NSentropy_[i] = Sn(temperature_[i]);
        }    

        if (progress != previousProgress) {
            postMessage(
            {messageType: "Progress", data: progress}
            );
            previousProgress = progress;
        }
    }   
  }



  


  function S(alpha, Tc, T){
    E1_ = 0;
    E2_ = 0.0015*T;
    dE_ = (E2_-E1_)/(NE_-1);
    let E = 0;
    let th1 = 0;
    let th2 = Math.PI/4.0;
    const Nth = 200;
    let dth = (th2-th1)/(Nth-1);
    let th = 0;
    let pg = 0;
    let f = 0;
    let fw = 0;
    let sum = 0;
    let del = 0;
    let d0 = delta(alpha,Tc,T);
        pg = Eg_*Math.cos(2.0*th1)
        del = d0*Math.cos(2.0*th1);    
        E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(pg,2));
        f = fermi(E,0,T);
        fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
        //if(f>0)
        if(!isNaN(fw))
        sum += 0.5*0.5*fw;

        E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(pg,2));
        f = fermi(E,0,T);
        //if(f>0)
        fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
        if(!isNaN(fw))
        sum += 0.5*0.5*fw;

        for(let i = 1; i<NE_-1; i++){  
            E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(pg,2));
            f = fermi(E,0,T);
            fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
            //if(f>0)
            if(!isNaN(fw))
            sum+=0.5*fw;
        }

        pg = Eg_*Math.cos(2.0*th2)
        del = d0*Math.cos(2.0*th2);    
        E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(pg,2));
        f = fermi(E,0,T);
        fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
        //if(f>0)
        if(!isNaN(fw))
        sum += 0.5*0.5*fw;

        E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(pg,2));
        f = fermi(E,0,T);
        fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
        //if(f>0)
        if(!isNaN(fw))
        sum += 0.5*0.5*fw;

        for(let i = 1; i<NE_-1; i++){  
            E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(pg,2));
            f = fermi(E,0,T);
            fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
            //if(f>0)
            if(!isNaN(fw))
            sum+=0.5*fw;
        }

    for(let j = 1; j<Nth-1; j++){
        th = th1+j*dth;
        pg = Eg_*Math.cos(2.0*th)
        del = d0*Math.cos(2.0*th);    
        E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(pg,2));
        f = fermi(E,0,T);
        fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
        //if(f>0)
        if(!isNaN(fw))
        sum += 0.5*fw;

        E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(pg,2));
        f = fermi(E,0,T);
        //if(f>0)
        fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
        if(!isNaN(fw))
        sum += 0.5*fw;

        for(let i = 1; i<NE_-1; i++){  
            E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(pg,2));
            f = fermi(E,0,T);
            fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f);
            if(!isNaN(fw))
            sum+=fw;
            //if(isNaN(sum)) console.log(f+' '+T);
        }
    }
     return sum*=(-dE_*dth*4.0/Math.PI*2.0*1E3*8.3145);

}



function Sn(T){
    return S(0,0,T);
}



function calculateSusceptibility(alpha,Tc){
	// Calculate the progress percentage.
    let progress = 0;//Math.round(i/list.length*100);  // Only send a progress update if the progress has changed
	
	for(let i = 0; i<NT_; i++){
		if(temperature_[i]==0){
            if(Eg_==0){
                NSsusceptibility_[i] = 5.788E-5*9.274E-24*4.0*Math.PI*1E-7*6.02E23;
            }else{
			    NSsusceptibility_[i] = 0;
            }
            if(alpha==0){
                susceptibility_[i] = NSsusceptibility_[i];
            }else{    
			    susceptibility_[i] = 0;
            }
		}else{
			susceptibility_[i] = Chi(alpha,Tc,temperature_[i]);
			NSsusceptibility_[i] = Chin(temperature_[i]);
		}    
		 if (progress != previousProgress) {
            postMessage(
            {messageType: "Progress", data: progress}
            );
            previousProgress = progress;
        }
	}
}


function Chi(alpha, Tc, T){
	E1_ = 0;
    E2_ = 0.0015*T;
    dE_ = (E2_-E1_)/(NE_-1);
    let E = 0;
    let th1 = 0;
    let th2 = Math.PI/4.0;
    const Nth = 200;
    let dth = (th2-th1)/(Nth-1);
    let th = 0;
    
    let df = 0;
    let sum = 0;
    let del = 0;
    let eg = 0;
    let d0 = delta(alpha,Tc,T);
    
        del = d0*Math.cos(2.0*th1);   
        eg = Eg_*Math.cos(2.0*th1); 
        //E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
        E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(eg,2));
        df = dfdE(E,T);
        if(!isNaN(df))
        sum += 0.5*0.5*df;


        //E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
        E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(eg,2));
        df = dfdE(E,T);
        if(!isNaN(df))
        sum += 0.5*0.5*df;

        for(let i = 1; i<NE_-1; i++){  
            //E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
            E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(eg,2));
            df = dfdE(E,T);
            if(!isNaN(df))
            sum+=0.5*df;
        }

        del = d0*Math.cos(2.0*th2);
        eg = Eg_*Math.cos(2.0*th2);      
        //E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
        E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(eg,2));
        df = dfdE(E,T);
        if(!isNaN(df))
        sum += 0.5*0.5*df;

        //E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
        E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(eg,2));
        df = dfdE(E,T);
        if(!isNaN(df))
        sum += 0.5*0.5*df;

        for(let i = 1; i<NE_-1; i++){  
            //E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
            E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(eg,2));
            df = dfdE(E,T);
            if(!isNaN(df))
            sum+=0.5*df;
        }

    for(let j = 1; j<Nth-1; j++){
        th = th1+j*dth;
        del = d0*Math.cos(2.0*th);    
        eg = Eg_*Math.cos(2.0*th);    
        //E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
        E = Math.sqrt(Math.pow(E1_,2)+Math.pow(del,2)+Math.pow(eg,2));
        df = dfdE(E,T);
        if(!isNaN(df))
        sum += 0.5*df;

        //E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
        E = Math.sqrt(Math.pow(E2_,2)+Math.pow(del,2)+Math.pow(eg,2));
        df = dfdE(E,T);
        if(!isNaN(df))
        sum += 0.5*df;

        for(let i = 1; i<NE_-1; i++){  
            //E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(Eg_*Math.cos(2.0*th),2));
            E = Math.sqrt(Math.pow(E1_+i*dE_,2)+Math.pow(del,2)+Math.pow(eg,2));
            df = dfdE(E,T);
            sum+=df;
            //if(isNaN(sum)) console.log(f+' '+T);
        }
    }  
		
	 return sum*=(-dE_*dth*4.0/Math.PI*2.0*5.788E-5*9.274E-24*4.0*Math.PI*1E-7*6.02E23); //muB = 5.788E-5 eV/T or 9.274E-24 J/T, mu0 = 4PiE-7, NA = 6.02E23. Should give susceptibility in m3/mol.

}

function Chin(T){
	return Chi(0,Tc_,T);
}


function calculateDelta(){
    for(let i=0; i<NT_; i++){
        gap_[i] = 1E3*delta(alpha_,Tc_,temperature_[i]);
    }
} 


function delta(alpha,Tc,T){
    let a = 7.0/5.0;
    let DCoverC = 0.955;
    let y = a*DCoverC*(Tc/T-1.0);
    if(y<0) return 0;
    return alpha*kb_*Tc*Math.tanh(Math.PI/alphaBCS_*Math.sqrt(y));
}


function fermi(E, mu, T){
    if(T==0){
        if(E<mu){
            return 1.0;
        }
        else if(E==mu){ 
            return 0.5;
        }
        else if(E>mu){
            return 0;
        }
    }
    return 1.0/(Math.exp((E-mu)/kb_/T)+1.0);
}


function dfdt(t,Delta,Eg,T){
    let x = Math.sqrt(Math.pow(t,2)+Math.pow(Delta,2)+Math.pow(Eg,2));
    let y = Math.exp(x/kb_/T);
    return -t*y/(Math.pow(y+1,2)*x*kb_*T);
}


function dfdE(E,T){
    let y = Math.exp(E/kb_/T);
    return -y/(Math.pow(y+1,2)*kb_*T);
}

function calculateLvsT(){
    let th1 = 0;
    let th2 = Math.PI/4.0;
    const Nth = 200;
    let dth = (th2-th1)/(Nth-1);
    let th = 0;
    let d0 = 0;
    for(let i = 0; i<NT_; i++){
        d0 = delta(alpha_,Tc_,temperature_[i]);
        for(let j = 0; j<Nth; j++){
            th = th1+j*dth;
            lambda_[i] += calculateLambda(0.02,d0*Math.cos(2.0*th),Eg_*Math.cos(2.0*th),temperature_[i]);
        }
        lambda_[i]*=dth/Math.PI*4.0;    
    }
}



function calculateLambda(mu, Delta, Eg, T){
    E1_ = -0.1;
    E2_ = -E1_;
    let NEL = 2*NE_;
    dE_ = (E2_-E1_)/(NEL-1);
    let E = 0;
    let a = 0;
    let b = 0;
    let c = 0;
    let d = 0;
    let x = 0;
    let g = Math.pow(Delta,2)+Math.pow(Eg,2);
    let sum = 0;

    E = E1_;
        if(E==0){                
            sum+=-0.5*2.0*mu*Math.pow(Delta,2)*Math.exp(Math.sqrt(g)/kb_/T)/Math.pow(Math.exp(Math.sqrt(g)/kb_/T)+1.0,2)/kb_/T/g+0.5*mu/Math.pow(Delta,2)/Math.pow(g,1.5)*(1.0-2.0*fermi(Math.sqrt(g),0,T));
        }else{
        x = Math.sqrt(Math.pow(E,2)+g);
        a = Math.pow(Delta,2)/x;
        b = mu*Math.pow(Delta,2)/E/x;
        c = 2.0*dfdt(E,Delta,Eg,T);
        if(isNaN(c)) c=0;
        d = E*(1.0-2.0*fermi(x,0,T))/Math.pow(x,2);
        sum+=0.5*(a+b)*(c+d);
        }

        E = E2_;
        if(E==0){
            sum+=-0.5*2.0*mu*Math.pow(Delta,2)*Math.exp(Math.sqrt(g)/kb_/T)/Math.pow(Math.exp(Math.sqrt(g)/kb_/T)+1.0,2)/kb_/T/g+0.5*mu/Math.pow(Delta,2)/Math.pow(g,1.5)*(1.0-2.0*fermi(Math.sqrt(g),0,T));
        }else{
        x = Math.sqrt(Math.pow(E,2)+g);
        a = Math.pow(Delta,2)/x;
        b = mu*Math.pow(Delta,2)/E/x;
        c = 2.0*dfdt(E,Delta,Eg,T);
        if(isNaN(c)) c=0;
        d = E*(1.0-2.0*fermi(x,0,T))/Math.pow(x,2);
        sum+=0.5*(a+b)*(c+d);
        }    

    for(let i = 1; i<NEL-1; i++){
        E = E1_+i*dE_;
        
        if(E==0){
            sum+=-2.0*mu*Math.pow(Delta,2)*Math.exp(Math.sqrt(g)/kb_/T)/Math.pow(Math.exp(Math.sqrt(g)/kb_/T)+1.0,2)/kb_/T/g+mu/Math.pow(Delta,2)/Math.pow(g,1.5)*(1.0-2.0*fermi(Math.sqrt(g),0,T));
        }else{
        x = Math.sqrt(Math.pow(E,2)+g);
        a = Math.pow(Delta,2)/x;
        b = mu*Math.pow(Delta,2)/E/x;
        c = 2.0*dfdt(E,Delta,Eg,T);
        if(isNaN(c)) c=0;
        d = E*(1.0-2.0*fermi(x,0,T))/Math.pow(x,2);
        if(isNaN(d)) console.log(d);
        sum+=(a+b)*(c+d);
        }
    }
    return sum*=dE_;
}
