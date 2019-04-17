#ifndef LINESEARCH_H
#define LINESEARCH_H
#include <Eigen/Core>

template<typename TProblem>
bool lineSearch(std::shared_ptr<TProblem> op, Eigen::VectorXd x, Eigen::MatrixXd &pos, Eigen::VectorXd dir, double &rate)
{
    Eigen::VectorXd fullX = op -> getFullVariables(x);
    int num = fullX.size() / 4;
    Eigen::VectorXd L = fullX.segment(0, 3 * num);
    Eigen::VectorXd S = fullX.segment(3*num, num);
                                      
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 1e-15;
    double maxStep = op->maxStep(x, dir);
    if(maxStep < alpha)
    {
        std::cout<<"Maximum line search step is too small: "<<maxStep<<std::endl;
        return false;
    }
    double beta = maxStep;

    std::cout<<"Maximum step is: "<<beta<<std::endl;
    
    if(rate > beta)
        rate = (alpha + beta) / 2;
    
    Eigen::VectorXd grad;
    double orig = op->value(L, S, pos);
    op->gradient(L, S, pos, grad);
    double deriv = dir.dot(grad);
    
    while (true)
    {
        Eigen::VectorXd newdE, newX, newFullX;
        newX =  x + rate*dir;
        newFullX = op -> getFullVariables(newX);
        Eigen::VectorXd newL = newFullX.segment(0, L.size());
        Eigen::MatrixXd newPos = pos;
        Eigen::VectorXd newS = newFullX.segment(L.size(), S.size());
        op->projectBack(newL, newS, newPos);
        double newenergy = op->value(newL, newS, newPos);
        op->gradient(newL, newS, newPos, newdE);
        Eigen::VectorXd fullGrad = op -> projM.transpose() * newdE;
        std::cout << "Trying rate = " << rate << ", energy now " << newenergy<<", L update "<<(newL-L).norm()<<", S update "<<(newS - S).norm() << ", L grad "<<fullGrad.segment(0, L.rows()).norm()<<" , S grad "<<fullGrad.segment(L.rows(), S.rows()).norm() << std::endl;
        
        if (std::isnan(newenergy) || newenergy > orig + rate*deriv*c1)
        {
            //            std::cout<<"Voilate the first Wolfe Condition"<<std::endl;
            beta = rate;
            rate = 0.5*(alpha + beta);
            if (beta - alpha < 1e-15)
            {
                rate = 1e-15;
                std::cout<<"Line Search failed, finished with Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                //                pos = newPos;
                return false;
            }
        }
        else if (newdE.dot(dir) < c2*deriv)
        {
            //            std::cout<<"Voilate the second Wolfe Condition"<<std::endl;
            alpha = rate;
            if (beta == maxStep)
            {
                rate = std::min(2 * alpha, beta);
            }
            else
            {
                rate = 0.5*(alpha + beta);
            }
            
            if (beta - alpha < 1e-10)
            {
                if(newenergy > orig)
                {
                    std::cout<<"Line Search failed with beta - alph < 1e-10 without decreasing energy, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                    return false;
                }
                else
                {
                    std::cout<<"Line Search succeed with beta - alph < 1e-10 without sufficient decrease, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                    return true;
                }
            }
        }
        else
        {
            std::cout<<"Line Search Finished with Rate = "<<rate<<std::endl;
            return true;
        }
    }
}


//template<typename TProblem>
//bool lineSearchForS(std::shared_ptr<TProblem> op, Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd &pos, Eigen::VectorXd dir, double &rate)
//{
//    double c1 = 0.1;
//    double c2 = 0.9;
//    double alpha = 1e-15;
//    double maxStep = op->maxStep(L,S,dir);
//    double beta = maxStep;
//
//    std::cout<<"Maximum step is: "<<beta<<std::endl;
//
//    if(rate > beta)
//        rate = (alpha + beta) / 2;
//
//    Eigen::VectorXd grad, gradS;
//    double orig = op->value(L, S, pos);
//    op->gradient(L, S, pos, grad);
//
//    gradS = grad.segment(L.size(), S.size());
//    double deriv = dir.dot(gradS);
//
//    while (true)
//    {
//        Eigen::VectorXd newdE, newdES;
//        Eigen::MatrixXd newPos = pos;
//        Eigen::VectorXd newS = S + rate*dir;
//        op->projectBack(L, newS, newPos);
//        double newenergy = op->value(L, newS, newPos);
//        op->gradient(L, newS, newPos, newdE);
//
//        newdES = newdE.segment(L.size(), S.size());
//
//
//        if (std::isnan(newenergy) || newenergy > orig + rate*deriv*c1)
//        {
//            //            std::cout<<"Voilate the first Wolfe Condition"<<std::endl;
//            beta = rate;
//            rate = 0.5*(alpha + beta);
//            if (beta - alpha < 1e-15)
//            {
//                rate = 1e-15;
//                std::cout<<"Line Search failed, finished with Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                //                pos = newPos;
//                return false;
//            }
//        }
//        else if (newdES.dot(dir) < c2*deriv)
//        {
//            alpha = rate;
//            if (beta == maxStep)
//            {
//                rate = 2 * alpha;
//            }
//            else
//            {
//                rate = 0.5*(alpha + beta);
//            }
//
//            if (beta - alpha < 1e-10)
//            {
//                if(newenergy > orig)
//                {
//                    std::cout<<"Line Search failed with beta - alph < 1e-10 without decreasing energy, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                    return false;
//                }
//                else
//                {
//                    std::cout<<"Line Search succeed with beta - alph < 1e-10 without sufficient decrease, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                    return false;
//                }
//            }
//        }
//        else
//        {
//            std::cout<<"Line Search Finished with Rate = "<<rate<<std::endl;
//            return true;
//        }
//    }
//}
//
//
//template<typename TProblem>
//bool lineSearchForL(std::shared_ptr<TProblem> op, Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd &pos, Eigen::VectorXd dir, double &rate)
//{
//    double c1 = 1e-4;
//    double c2 = 0.9;
//    double alpha = 1e-15;
//    double maxStep = 1e15;
//    double beta = maxStep;
//
//    std::cout<<"Maximum step is: "<<beta<<std::endl;
//
//    if(rate > beta)
//        rate = (alpha + beta) / 2;
//
//    Eigen::VectorXd grad, gradL;
//    double orig = op->value(L, S, pos);
//    op->gradient(L, S, pos, grad);
//
//    gradL = grad.segment(0, L.size());
//    double deriv = dir.dot(gradL);
//
//    while (true)
//    {
//        Eigen::VectorXd newdE, newdEL;
//        Eigen::MatrixXd newPos = pos;
//        Eigen::VectorXd newL = L + rate*dir;
//        op->projectBack(newL, S, newPos);
//        double newenergy = op->value(newL, S, newPos);
//        op->gradient(newL, S, newPos, newdE);
//
//        newdEL = newdE.segment(0, L.size());
//
//
//        if (std::isnan(newenergy) || newenergy > orig + rate*deriv*c1)
//        {
//            //            std::cout<<"Voilate the first Wolfe Condition"<<std::endl;
//            beta = rate;
//            rate = 0.5*(alpha + beta);
//            if (beta - alpha < 1e-15)
//            {
//                rate = 1e-15;
//                std::cout<<"Line Search failed, finished with Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                //                pos = newPos;
//                return false;
//            }
//        }
//        else if (newdEL.dot(dir) < c2*deriv)
//        {
//            alpha = rate;
//            if (beta == maxStep)
//            {
//                rate = 2 * alpha;
//            }
//            else
//            {
//                rate = 0.5*(alpha + beta);
//            }
//
//            if (beta - alpha < 1e-10)
//            {
//                if(newenergy > orig)
//                {
//                    std::cout<<"Line Search failed with beta - alph < 1e-10 without decreasing energy, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                    return false;
//                }
//                else
//                {
//                    std::cout<<"Line Search succeed with beta - alph < 1e-10 without sufficient decrease, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
//                    return false;
//                }
//            }
//        }
//        else
//        {
//            std::cout<<"Line Search Finished with Rate = "<<rate<<std::endl;
//            return true;
//        }
//    }
//}

#endif
