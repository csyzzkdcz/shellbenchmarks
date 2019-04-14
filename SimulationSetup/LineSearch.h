#include <Eigen/Core>

template<typename TProblem>
bool lineSearch(std::shared_ptr<TProblem> op, Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd &pos, Eigen::VectorXd dir, double &rate)
{
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 1e-15;
    double maxStep = op->maxStep(L,S,dir);
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
        Eigen::VectorXd newdE;
        Eigen::VectorXd newL = L + rate*dir.segment(0, L.size());
        Eigen::MatrixXd newPos = pos;
        Eigen::VectorXd newS = S + rate*dir.segment(L.size(), S.size());
        op->projectBack(newL, newS, newPos);
        double newenergy = op->value(newL, newS, newPos);
        op->gradient(newL, newS, newPos, newdE);
        
        std::cout << "Trying rate = " << rate << ", energy now " << newenergy<<", L update "<<(newL-L).norm()<<", S update "<<(newS - S).norm() << ", L grad "<<newdE.segment(0, L.rows()).norm()<<" , S grad "<<newdE.segment(L.rows(), S.rows()).norm() << std::endl;
        
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
                rate = 2 * alpha;
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
                    return false;
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
